/* Copyright (C) 2005-2025 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* This file implements dispersive materials for Meep via a
   polarization P = \chi(\omega) W, where W is e.g. E or H.  Each
   subclass of the susceptibility class should implement a different
   type of \chi(\omega).  The subclass knows how to timestep P given W
   at the current (and possibly previous) timestep, and any additional
   internal data that needs to be allocated along with P.

   Each \chi(\omega) is spatially multiplied by a (scalar) sigma
   array.  The meep::fields class is responsible for allocating P and
   sigma and passing them to susceptibility::update_P. */

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <numeric>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netdb.h>
#include <unordered_map>
#include <unistd.h>
#include "meep.hpp"
#include "meep_internals.hpp"
#include <iostream>
#include <chrono>
#include <support/ziggurat.hpp>
#include <support/pcg_random.hpp>
#include <support/pcg_extras.hpp>
#include <support/pcg_uint128.hpp>

using namespace std;
using namespace cxx;

namespace meep {

int susceptibility::cur_id = 0;

susceptibility *susceptibility::clone() const {
  susceptibility *sus = new susceptibility(*this);
  sus->next = 0;
  sus->ntot = ntot;
  sus->id = id;
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) {
    if (sigma[c][d]) {
      sus->sigma[c][d] = new realnum[ntot];
      memcpy(sus->sigma[c][d], sigma[c][d], sizeof(realnum) * ntot);
    }
    else
      sus->sigma[c][d] = NULL;
    sus->trivial_sigma[c][d] = trivial_sigma[c][d];
  }
  return sus;
}

// generic base class definition.
std::complex<realnum> susceptibility::chi1(realnum freq, realnum sigma) {
  (void)freq;
  (void)sigma;
  return std::complex<realnum>(0, 0);
}

void susceptibility::delete_internal_data(void *data) const { free(data); }

/* Return whether or not we need to allocate P[c][cmp].  (We don't need to
   allocate P[c] if we can be sure it will be zero.)

   We are a bit wasteful because if sigma is nontrivial in *any* chunk,
   we allocate the corresponding P on *every* owned chunk.  This greatly
   simplifies communication in boundaries.cpp, because we can be sure that
   one chunk has a P then any chunk it borders has the same P, so we don't
   have to worry about communicating with something that doesn't exist.
   TODO: reduce memory usage (bookkeeping seem much harder, though).
*/
bool susceptibility::needs_P(component c, int cmp, realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  if (!is_electric(c) && !is_magnetic(c)) return false;
  FOR_DIRECTIONS(d) {
    if (!trivial_sigma[c][d] && W[direction_component(c, d)][cmp]) return true;
  }
  return false;
}

/* return whether we need the notowned parts of the W field --
   by default, this is only the case if sigma has offdiagonal components
   coupling P to W.   (See needs_P: again, this true if the notowned
   W is needed in *any* chunk.) */
bool susceptibility::needs_W_notowned(component c, realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  FOR_DIRECTIONS(d) {
    if (d != component_direction(c)) {
      component cP = direction_component(c, d);
      if (needs_P(cP, 0, W) && !trivial_sigma[cP][component_direction(c)]) return true;
    }
  }
  return false;
}

typedef struct {
  size_t sz_data;
  size_t ntot;
  realnum *P[NUM_FIELD_COMPONENTS][2];
  realnum *P_prev[NUM_FIELD_COMPONENTS][2];
  realnum data[1];
} lorentzian_data;

// for Lorentzian susc. the internal data is just a backup of P from
// the previous timestep.
void *lorentzian_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
                                                   const grid_volume &gv) const {
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) num += 2 * gv.ntot();
  }
  size_t sz = sizeof(lorentzian_data) + sizeof(realnum) * (num - 1);
  lorentzian_data *d = (lorentzian_data *)malloc(sz);
  if (d == NULL) meep::abort("%s:%i:out of memory(%lu)", __FILE__, __LINE__, sz);
  d->sz_data = sz;
  return (void *)d;
}

void lorentzian_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], realnum dt,
                                                   const grid_volume &gv, void *data) const {
  (void)dt; // unused
  lorentzian_data *d = (lorentzian_data *)data;
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  size_t ntot = d->ntot = gv.ntot();
  realnum *P = d->data;
  realnum *P_prev = d->data + ntot;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) {
      d->P[c][cmp] = P;
      d->P_prev[c][cmp] = P_prev;
      P += 2 * ntot;
      P_prev += 2 * ntot;
    }
  }
}

void *lorentzian_susceptibility::copy_internal_data(void *data) const {
  lorentzian_data *d = (lorentzian_data *)data;
  if (!d) return 0;
  lorentzian_data *dnew = (lorentzian_data *)malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);
  size_t ntot = d->ntot;
  realnum *P = dnew->data;
  realnum *P_prev = dnew->data + ntot;
  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      dnew->P[c][cmp] = P;
      dnew->P_prev[c][cmp] = P_prev;
      P += 2 * ntot;
      P_prev += 2 * ntot;
    }
  }
  return (void *)dnew;
}

#if 0
/* Return true if the discretized Lorentzian ODE is intrinsically unstable,
   i.e. if it corresponds to a filter with a pole z outside the unit circle.
   Note that the pole satisfies the quadratic equation:
            (z + 1/z - 2)/dt^2 + g*(z - 1/z)/(2*dt) + w^2 = 0
   where w = 2*pi*omega_0 and g = 2*pi*gamma.   It is just a little
   algebra from this to get the condition for a root with |z| > 1.

   FIXME: this test seems to be too conservative (issue #12) */
static bool lorentzian_unstable(realnum omega_0, realnum gamma, realnum dt) {
  realnum w = 2 * pi * omega_0, g = 2 * pi * gamma;
  realnum g2 = g * dt / 2, w2 = (w * dt) * (w * dt);
  realnum b = (1 - w2 / 2) / (1 + g2), c = (1 - g2) / (1 + g2);
  return b * b > c && 2 * b * b - c + 2 * fabs(b) * sqrt(b * b - c) > 1;
}
#endif

#define SWAP(t, a, b)                                                                              \
  {                                                                                                \
    t SWAP_temp = a;                                                                               \
    a = b;                                                                                         \
    b = SWAP_temp;                                                                                 \
  }

// stable averaging of offdiagonal components
#define OFFDIAG(u, g, sx, s)                                                                       \
  (0.25 * ((g[i] + g[i - sx]) * u[i] + (g[i + s] + g[(i + s) - sx]) * u[i + s]))

namespace {

const double mxl_fs_to_au = 41.341373335;
const double mxl_efield_mu_to_au_prefactor = 1.2929541569381223e-6;
const double mxl_source_amp_au_to_mu = 0.002209799779149953;
const int mxl_molecules_per_chunk_id_block = 100000;
const int mxl_chunks_per_rank_id_block = 100;
const int mxl_max_molecule_id = 2147483647;
const char *mxl_socket_mapping_filename = "mxl_socket_molecule_map.csv";

/* These counters are reset whenever material susceptibilities are rebuilt. */
int mxl_next_chunk_ordinal = 0;
size_t mxl_local_active_site_count = 0;
size_t mxl_local_required_driver_count = 0;
size_t mxl_expected_total_required_driver_count = 0;
bool mxl_socket_susceptibility_present = false;
bool mxl_driver_count_reported = false;
bool mxl_socket_imag_field_coupling_active = false;

struct mxl_socket_mapping_record {
  int molecule_id;
  int chunk_ordinal;
  int drive_cmp;
  size_t site_ordinal;
  size_t local_index;
  double x_or_r;
  double y;
  double z;
  std::string label;
  std::string host;
  int port;
};

std::vector<mxl_socket_mapping_record> mxl_local_mapping_records;

static int mxl_dim_power(ndim dim) {
  switch (dim) {
    case D1: return 1;
    case D2: return 2;
    case D3: return 3;
    /* The local polarization-density update remains normalized by the
       equivalent 3d Cartesian Yee-cell volume.  Cylindrical annular
       integration weights belong in global integrals, not in the local
       polarization-density update; rescaling_factor is a symmetric
       bright-state coupling scale applied to drive and response. */
    case Dcyl: return 3;
  }
  return 3;
}

static void mxl_field_components(ndim dim, component comps[3]) {
  if (dim == Dcyl) {
    /* The socket protocol remains Cartesian.  Cylindrical Meep fields are
       modal coefficients in the local (r,phi,z) basis; at phi=0 this maps
       to the socket's (x,y,z) slots. */
    comps[0] = Er;
    comps[1] = Ep;
    comps[2] = Ez;
  }
  else {
    comps[0] = Ex;
    comps[1] = Ey;
    comps[2] = Ez;
  }
}

static int mxl_amp_axis(component c) {
  switch (component_direction(c)) {
    case X:
    case R: return 0;
    case Y:
    case P: return 1;
    case Z: return 2;
    default: return 0;
  }
}

static void mxl_assert_little_endian() {
  const uint16_t x = 1;
  if (*((const unsigned char *)&x) != 1)
    meep::abort("MXLSocketSusceptibility requires a little-endian host.");
}

static void mxl_append_i32(std::vector<unsigned char> &buf, int value) {
  uint32_t v = (uint32_t)value;
  buf.push_back((unsigned char)(v & 0xffu));
  buf.push_back((unsigned char)((v >> 8) & 0xffu));
  buf.push_back((unsigned char)((v >> 16) & 0xffu));
  buf.push_back((unsigned char)((v >> 24) & 0xffu));
}

static int mxl_read_i32(const unsigned char *p) {
  uint32_t v = ((uint32_t)p[0]) | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) |
               ((uint32_t)p[3] << 24);
  return (int)v;
}

static void mxl_write_f64(unsigned char *p, double value) { memcpy(p, &value, sizeof(double)); }

static double mxl_read_f64(const unsigned char *p) {
  double value;
  memcpy(&value, p, sizeof(double));
  return value;
}

static double mxl_x_or_r(const vec &r) {
  if (r.dim == Dcyl) return r.in_direction(R);
  return has_direction(r.dim, X) ? r.in_direction(X) : 0.0;
}

static double mxl_y(const vec &r) { return has_direction(r.dim, Y) ? r.in_direction(Y) : 0.0; }

static double mxl_z(const vec &r) { return has_direction(r.dim, Z) ? r.in_direction(Z) : 0.0; }

static void mxl_fprint_csv_string(FILE *f, const std::string &s) {
  fputc('"', f);
  for (size_t i = 0; i < s.size(); ++i) {
    if (s[i] == '"') fputc('"', f);
    fputc(s[i], f);
  }
  fputc('"', f);
}

static void mxl_write_mapping_csv() {
  begin_critical_section(31416);

  FILE *f = fopen(mxl_socket_mapping_filename, my_rank() == 0 ? "w" : "a");
  if (!f)
    meep::abort("MXLSocketSusceptibility failed to open molecule map '%s' for writing: %s",
                mxl_socket_mapping_filename, strerror(errno));

  if (my_rank() == 0) {
    fprintf(f, "molecule_id,rank,chunk_ordinal,label,host,port,drive_cmp,drive_part,"
               "site_ordinal,local_index,x_or_r,y,z\n");
  }

  for (size_t i = 0; i < mxl_local_mapping_records.size(); ++i) {
    const mxl_socket_mapping_record &r = mxl_local_mapping_records[i];
    fprintf(f, "%d,%d,%d,", r.molecule_id, my_rank(), r.chunk_ordinal);
    mxl_fprint_csv_string(f, r.label);
    fprintf(f, ",");
    mxl_fprint_csv_string(f, r.host);
    fprintf(f, ",%d,%d,%s,%zu,%zu,%.17g,%.17g,%.17g\n", r.port, r.drive_cmp,
            r.drive_cmp == 0 ? "real" : "imag", r.site_ordinal, r.local_index, r.x_or_r, r.y,
            r.z);
  }

  if (fclose(f) != 0)
    meep::abort("MXLSocketSusceptibility failed to close molecule map '%s': %s",
                mxl_socket_mapping_filename, strerror(errno));

  end_critical_section(31416);
  all_wait();

  master_printf("MXLSocketSusceptibility molecule map written to %s\n",
                mxl_socket_mapping_filename);
}

/* Minimal TCP client for the aggregate MaxwellLink susceptibility protocol. */
class mxl_socket_client {
public:
  mxl_socket_client() : fd(-1), cached_nreq(0), cached_molecule_id_data(NULL) {}
  mxl_socket_client(const mxl_socket_client &)
      : fd(-1), cached_nreq(0), cached_molecule_id_data(NULL) {}
  mxl_socket_client &operator=(const mxl_socket_client &) {
    close();
    clear_step_cache();
    fd = -1;
    return *this;
  }
  ~mxl_socket_client() { close(); }

  void ensure_connected(const std::string &host, int port, double timeout) {
    if (fd >= 0) return;

    std::ostringstream port_stream;
    port_stream << port;

    struct addrinfo hints;
    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;

    struct addrinfo *result = NULL;
    int rc = getaddrinfo(host.c_str(), port_stream.str().c_str(), &hints, &result);
    if (rc != 0)
      meep::abort("MXLSocketSusceptibility getaddrinfo(%s:%d) failed: %s", host.c_str(), port,
                  gai_strerror(rc));

    int connected_fd = -1;
    for (struct addrinfo *rp = result; rp != NULL; rp = rp->ai_next) {
      connected_fd = socket(rp->ai_family, rp->ai_socktype, rp->ai_protocol);
      if (connected_fd < 0) continue;
      set_timeout(connected_fd, timeout);
      if (connect(connected_fd, rp->ai_addr, rp->ai_addrlen) == 0) break;
      ::close(connected_fd);
      connected_fd = -1;
    }

    freeaddrinfo(result);

    if (connected_fd < 0)
      meep::abort("MXLSocketSusceptibility failed to connect to MaxwellLink endpoint %s:%d",
                  host.c_str(), port);
    fd = connected_fd;
  }

  void send_init(const std::vector<int> &molecule_ids, double dt_au, double rescaling_factor,
                 double time_units_fs, double timeout, int rank, size_t expected_total_molecules) {
    std::ostringstream json;
    json.precision(17);
    json << "{\"protocol\":\"mxl_socket_susceptibility_v1\",";
    json << "\"dt_au\":" << dt_au << ",";
    json << "\"rank\":" << rank << ",";
    json << "\"rescaling_factor\":" << rescaling_factor << ",";
    json << "\"time_units_fs\":" << time_units_fs << ",";
    json << "\"timeout\":" << timeout << ",";
    json << "\"expected_total_molecules\":" << expected_total_molecules << ",";
    json << "\"molecule_ids\":[";
    for (size_t i = 0; i < molecule_ids.size(); ++i) {
      if (i) json << ",";
      json << molecule_ids[i];
    }
    json << "]}";

    send_msg("MXLINIT");
    send_bytes(json.str());
    std::string ready = recv_msg();
    if (ready != "MXLREADY")
      meep::abort("MXLSocketSusceptibility expected MXLREADY after MXLINIT, got '%s'",
                  ready.c_str());
  }

  void step(const std::vector<int> &molecule_ids, const std::vector<double> &efields_au,
            std::vector<double> &amps_au) {
    const size_t nreq = molecule_ids.size();
    if (efields_au.size() != 3 * nreq)
      meep::abort("MXLSocketSusceptibility internal efield array size mismatch.");

    ensure_step_cache(molecule_ids);
    unsigned char *field_data = step_frame.data() + 20;
    for (size_t i = 0; i < nreq; ++i) {
      mxl_write_f64(field_data + 24 * i + 0, efields_au[3 * i + 0]);
      mxl_write_f64(field_data + 24 * i + 8, efields_au[3 * i + 1]);
      mxl_write_f64(field_data + 24 * i + 16, efields_au[3 * i + 2]);
    }
    send_all(step_frame.data(), step_frame.size());

    recv_all(result_head.data(), result_head.size());
    std::string header = parse_header(result_head.data());
    if (header != "AGGRESULT")
      meep::abort("MXLSocketSusceptibility expected AGGRESULT, got '%s'", header.c_str());

    int nresp = mxl_read_i32(result_head.data() + 12);
    if (nresp < 0) meep::abort("MXLSocketSusceptibility received negative response count.");

    result_fixed.resize((size_t)nresp * 32);
    if (!result_fixed.empty()) recv_all(result_fixed.data(), result_fixed.size());

    if (amps_au.size() != 3 * nreq) amps_au.resize(3 * nreq);
    std::fill(amps_au.begin(), amps_au.end(), 0.0);
    std::fill(seen.begin(), seen.end(), 0);
    size_t total_extra = 0;

    for (int j = 0; j < nresp; ++j) {
      const unsigned char *rec = result_fixed.data() + (size_t)j * 32;
      int mid = mxl_read_i32(rec);
      size_t idx;
      if ((size_t)j < nreq && molecule_ids[(size_t)j] == mid) {
        idx = (size_t)j;
      }
      else {
        std::unordered_map<int, size_t>::const_iterator it = id_to_index.find(mid);
        if (it == id_to_index.end())
          meep::abort("MXLSocketSusceptibility received response for unknown molecule_id %d", mid);
        idx = it->second;
      }
      amps_au[3 * idx + 0] = mxl_read_f64(rec + 4);
      amps_au[3 * idx + 1] = mxl_read_f64(rec + 12);
      amps_au[3 * idx + 2] = mxl_read_f64(rec + 20);
      int extra_len = mxl_read_i32(rec + 28);
      if (extra_len < 0)
        meep::abort("MXLSocketSusceptibility received negative extra payload length.");
      total_extra += (size_t)extra_len;
      seen[idx] = 1;
    }

    if (total_extra) {
      if (result_extras.size() < total_extra) result_extras.resize(total_extra);
      recv_all(result_extras.data(), total_extra);
    }
    for (size_t i = 0; i < nreq; ++i)
      if (!seen[i])
        meep::abort("MXLSocketSusceptibility missing response for molecule_id %d",
                    molecule_ids[i]);
  }

private:
  int fd;
  size_t cached_nreq;
  const int *cached_molecule_id_data;
  std::vector<int> cached_molecule_ids;
  std::unordered_map<int, size_t> id_to_index;
  std::vector<unsigned char> step_frame;
  std::vector<unsigned char> result_head;
  std::vector<unsigned char> result_fixed;
  std::vector<unsigned char> result_extras;
  std::vector<char> seen;

  static void set_timeout(int sockfd, double timeout) {
    if (timeout <= 0) return;
    struct timeval tv;
    tv.tv_sec = (time_t)timeout;
    tv.tv_usec = (suseconds_t)((timeout - tv.tv_sec) * 1000000.0);
    setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    setsockopt(sockfd, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
  }

  void close() {
    if (fd >= 0) {
      ::close(fd);
      fd = -1;
    }
  }

  void clear_step_cache() {
    cached_nreq = 0;
    cached_molecule_id_data = NULL;
    cached_molecule_ids.clear();
    id_to_index.clear();
    step_frame.clear();
    result_head.clear();
    result_fixed.clear();
    result_extras.clear();
    seen.clear();
  }

  void ensure_step_cache(const std::vector<int> &molecule_ids) {
    const int *molecule_id_data = molecule_ids.empty() ? NULL : &molecule_ids[0];
    /* The caller builds molecule_ids during susceptibility initialization and
       then reuses the same vector for every time step of this client. */
    if (!step_frame.empty() && cached_molecule_id_data == molecule_id_data &&
        cached_nreq == molecule_ids.size())
      return;

    cached_molecule_ids = molecule_ids;
    cached_molecule_id_data = molecule_id_data;
    cached_nreq = cached_molecule_ids.size();
    if (cached_nreq > (size_t)mxl_max_molecule_id)
      meep::abort("MXLSocketSusceptibility request count exceeds int32 protocol range.");

    step_frame.clear();
    step_frame.reserve(20 + 32 * cached_nreq);
    append_header(step_frame, "AGGSTEP");
    mxl_append_i32(step_frame, (int)cached_nreq);
    mxl_append_i32(step_frame, (int)cached_nreq);
    step_frame.resize(step_frame.size() + 24 * cached_nreq, 0);
    for (size_t i = 0; i < cached_nreq; ++i) {
      mxl_append_i32(step_frame, cached_molecule_ids[i]);
      mxl_append_i32(step_frame, (int)i);
    }

    id_to_index.clear();
    id_to_index.reserve(cached_nreq);
    for (size_t i = 0; i < cached_nreq; ++i)
      id_to_index[cached_molecule_ids[i]] = i;

    result_head.assign(16, 0);
    result_fixed.clear();
    seen.assign(cached_nreq, 0);
  }

  static void append_header(std::vector<unsigned char> &buf, const char *msg) {
    size_t len = strlen(msg);
    if (len > 12) meep::abort("MXLSocketSusceptibility internal header too long.");
    for (size_t i = 0; i < 12; ++i)
      buf.push_back(i < len ? (unsigned char)msg[i] : (unsigned char)' ');
  }

  static std::string parse_header(const unsigned char *p) {
    std::string s((const char *)p, 12);
    while (!s.empty() && s[s.size() - 1] == ' ')
      s.erase(s.size() - 1);
    return s;
  }

  void send_msg(const char *msg) {
    unsigned char hdr[12];
    memset(hdr, ' ', sizeof(hdr));
    size_t len = strlen(msg);
    if (len > sizeof(hdr)) meep::abort("MXLSocketSusceptibility internal header too long.");
    memcpy(hdr, msg, len);
    send_all(hdr, sizeof(hdr));
  }

  std::string recv_msg() {
    unsigned char hdr[12];
    recv_all(hdr, sizeof(hdr));
    return parse_header(hdr);
  }

  void send_bytes(const std::string &payload) {
    std::vector<unsigned char> lenbuf;
    mxl_append_i32(lenbuf, (int)payload.size());
    send_all(lenbuf.data(), lenbuf.size());
    if (!payload.empty()) send_all((const unsigned char *)payload.data(), payload.size());
  }

  void send_all(const unsigned char *buf, size_t len) {
    size_t done = 0;
    while (done < len) {
      ssize_t n = ::send(fd, buf + done, len - done, 0);
      if (n <= 0)
        meep::abort("MXLSocketSusceptibility socket send failed: %s", strerror(errno));
      done += (size_t)n;
    }
  }

  void recv_all(unsigned char *buf, size_t len) {
    size_t done = 0;
    while (done < len) {
      ssize_t n = ::recv(fd, buf + done, len - done, 0);
      if (n <= 0)
        meep::abort("MXLSocketSusceptibility socket receive failed: %s", strerror(errno));
      done += (size_t)n;
    }
  }
};

/* Per-chunk state for the socket-backed polarization variables. */
struct mxl_socket_data {
  size_t ntot;
  realnum *P[NUM_FIELD_COMPONENTS][2];
  std::vector<realnum> P_store[NUM_FIELD_COMPONENTS][2];
  std::vector<size_t> active_indices;
  std::vector<unsigned char> comp_mask; // owned nonzero-sigma bits (x/r=1, y/p=2, z=4)
  std::vector<int> drive_cmps;
  std::vector<int> molecule_ids;
  std::vector<double> efields_au;
  std::vector<double> amps_au;
  bool initialized;
  mxl_socket_client client;

  mxl_socket_data() : ntot(0), initialized(false), client() {
    FOR_COMPONENTS(c) DOCMP2 { P[c][cmp] = NULL; }
  }

  mxl_socket_data(const mxl_socket_data &other)
      : ntot(other.ntot), active_indices(other.active_indices), comp_mask(other.comp_mask),
        drive_cmps(other.drive_cmps),
        molecule_ids(other.molecule_ids), efields_au(other.efields_au),
        amps_au(other.amps_au), initialized(false), client() {
    FOR_COMPONENTS(c) DOCMP2 { P_store[c][cmp] = other.P_store[c][cmp]; }
    reset_ptrs();
  }

  void reset_ptrs() {
    FOR_COMPONENTS(c) DOCMP2 {
      P[c][cmp] = P_store[c][cmp].empty() ? NULL : P_store[c][cmp].data();
    }
  }
};

} // namespace

void reset_mxl_socket_susceptibility_chunk_ordinals() {
  mxl_next_chunk_ordinal = 0;
  mxl_local_active_site_count = 0;
  mxl_local_required_driver_count = 0;
  mxl_expected_total_required_driver_count = 0;
  mxl_local_mapping_records.clear();
  mxl_socket_susceptibility_present = false;
  mxl_driver_count_reported = false;
  mxl_socket_imag_field_coupling_active = false;
}

void report_mxl_socket_susceptibility_driver_count(field_type ft) {
  if (ft != E_stuff || mxl_driver_count_reported) return;

  /* This is called after all owned chunks have initialized their susceptibility
     data, but before any socket connection is opened. */
  const bool present = or_to_all(mxl_socket_susceptibility_present);
  if (!present) {
    mxl_driver_count_reported = true;
    return;
  }

  const size_t total_active_site_count = sum_to_all(mxl_local_active_site_count);
  const size_t total_required_driver_count = sum_to_all(mxl_local_required_driver_count);
  mxl_expected_total_required_driver_count = total_required_driver_count;
  const bool imag_field_coupling_active = or_to_all(mxl_socket_imag_field_coupling_active);
  begin_critical_section(31415);
  printf("MXLSocketSusceptibility rank %d: required socket drivers = %zu\n", my_rank(),
         mxl_local_required_driver_count);
  fflush(stdout);
  end_critical_section(31415);
  all_wait();
  master_printf("MXLSocketSusceptibility total required socket drivers = %zu\n",
                total_required_driver_count);
  master_printf("MXLSocketSusceptibility total active grid points = %zu\n",
                total_active_site_count);
  if (imag_field_coupling_active)
    master_printf("MXLSocketSusceptibility complex-field imaginary coupling enabled: "
                  "independent socket molecules are used for real and imaginary E-field "
                  "components, so required socket drivers are doubled relative to active "
                  "grid points.\n");
  if (total_required_driver_count > 0) mxl_write_mapping_csv();
  mxl_driver_count_reported = true;
}

void lorentzian_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], realnum dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *)P_internal_data;
  const realnum omega2pi = 2 * pi * omega_0, g2pi = gamma * 2 * pi;
  const realnum omega0dtsqr = omega2pi * omega2pi * dt * dt;
  const realnum gamma1inv = 1 / (1 + g2pi * dt / 2), gamma1 = (1 - g2pi * dt / 2);
  const realnum omega0dtsqr_denom = no_omega_0_denominator ? 0 : omega0dtsqr;
  (void)W_prev; // unused;

  // TODO: add back lorentzian_unstable(omega_0, gamma, dt) if we can improve the stability test
  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
      if (w && s) {
        realnum *p = d->P[c][cmp], *pp = d->P_prev[c][cmp];

        // directions/strides for offdiagonal terms, similar to update_eh
        const direction d = component_direction(c);
        const ptrdiff_t is = gv.stride(d) * (is_magnetic(c) ? -1 : +1);
        direction d1 = cycle_direction(gv.dim, d, 1);
        component c1 = direction_component(c, d1);
        ptrdiff_t is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
        const realnum *w1 = W[c1][cmp];
        const realnum *s1 = w1 ? sigma[c][d1] : NULL;
        direction d2 = cycle_direction(gv.dim, d, 2);
        component c2 = direction_component(c, d2);
        ptrdiff_t is2 = gv.stride(d2) * (is_magnetic(c) ? -1 : +1);
        const realnum *w2 = W[c2][cmp];
        const realnum *s2 = w2 ? sigma[c][d2] : NULL;

        if (s2 && !s1) { // make s1 the non-NULL one if possible
          SWAP(direction, d1, d2);
          SWAP(component, c1, c2);
          SWAP(ptrdiff_t, is1, is2);
          SWAP(const realnum *, w1, w2);
          SWAP(const realnum *, s1, s2);
        }
        if (s1 && s2) { // 3x3 anisotropic
          PLOOP_OVER_VOL_OWNED(gv, c, i) {
            // s[i] != 0 check is a bit of a hack to work around
            // some instabilities that occur near the boundaries
            // of materials; see PR #666
            if (s[i] != 0) {
              realnum pcur = p[i];
              p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] +
                                  omega0dtsqr * (s[i] * w[i] + OFFDIAG(s1, w1, is1, is) +
                                                 OFFDIAG(s2, w2, is2, is)));
              pp[i] = pcur;
            }
          }
        }
        else if (s1) { // 2x2 anisotropic
          PLOOP_OVER_VOL_OWNED(gv, c, i) {
            if (s[i] != 0) { // see above
              realnum pcur = p[i];
              p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] +
                                  omega0dtsqr * (s[i] * w[i] + OFFDIAG(s1, w1, is1, is)));
              pp[i] = pcur;
            }
          }
        }
        else { // isotropic
          PLOOP_OVER_VOL_OWNED(gv, c, i) {
            realnum pcur = p[i];
            p[i] = gamma1inv *
                   (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] + omega0dtsqr * (s[i] * w[i]));
            pp[i] = pcur;
          }
        }
      }
    }
  }
}

void lorentzian_susceptibility::subtract_P(field_type ft,
                                           realnum *f_minus_p[NUM_FIELD_COMPONENTS][2],
                                           void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *)P_internal_data;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  size_t ntot = d->ntot;
  FOR_FT_COMPONENTS(ft, ec) DOCMP2 {
    if (d->P[ec][cmp]) {
      component dc = field_type_component(ft2, ec);
      if (f_minus_p[dc][cmp]) {
        realnum *p = d->P[ec][cmp];
        realnum *fmp = f_minus_p[dc][cmp];
        for (size_t i = 0; i < ntot; ++i)
          fmp[i] -= p[i];
      }
    }
  }
}

int lorentzian_susceptibility::num_cinternal_notowned_needed(component c,
                                                             void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *)P_internal_data;
  return d->P[c][0] ? 1 : 0;
}

realnum *lorentzian_susceptibility::cinternal_notowned_ptr(int inotowned, component c, int cmp,
                                                           int n, void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *)P_internal_data;
  (void)inotowned; // always = 0
  if (!d || !d->P[c][cmp]) return NULL;
  return d->P[c][cmp] + n;
}

std::complex<realnum> lorentzian_susceptibility::chi1(realnum freq, realnum sigma) {
  if (no_omega_0_denominator) {
    // Drude model
    return sigma * omega_0 * omega_0 / std::complex<realnum>(-freq * freq, -gamma * freq);
  }
  else {
    // Standard Lorentzian model
    return sigma * omega_0 * omega_0 /
           std::complex<realnum>(omega_0 * omega_0 - freq * freq, -gamma * freq);
  }
}

void lorentzian_susceptibility::dump_params(h5file *h5f, size_t *start) {
  size_t num_params = 5;
  size_t params_dims[1] = {num_params};
  realnum params_data[] = {4, (realnum)get_id(), omega_0, gamma, (realnum)no_omega_0_denominator};
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;
}

void noisy_lorentzian_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                               realnum *W_prev[NUM_FIELD_COMPONENTS][2], realnum dt,
                                               const grid_volume &gv, void *P_internal_data) const {
  lorentzian_susceptibility::update_P(W, W_prev, dt, gv, P_internal_data);
  lorentzian_data *d = (lorentzian_data *)P_internal_data;

  const realnum g2pi = gamma * 2 * pi;
  const realnum w2pi = omega_0 * 2 * pi;
  const realnum amp = w2pi * noise_amp * sqrt(g2pi) * dt * dt / (1 + g2pi * dt / 2);
  /* for uniform random numbers in [-amp,amp] below, multiply amp by sqrt(3) */

  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      const realnum *s = sigma[c][component_direction(c)];
      if (s) {
        realnum *p = d->P[c][cmp];
        LOOP_OVER_VOL_OWNED(gv, c, i) { p[i] += gaussian_random(0, amp * sqrt(s[i])); }
        // for uniform random numbers, use uniform_random(-1,1) * amp * sqrt(s[i])
        // for gaussian random numbers, use gaussian_random(0, amp * sqrt(s[i]))
      }
    }
  }
}

void noisy_lorentzian_susceptibility::dump_params(h5file *h5f, size_t *start) {
  size_t num_params = 6;
  size_t params_dims[1] = {num_params};
  realnum params_data[] = {
      5, (realnum)get_id(), noise_amp, omega_0, gamma, (realnum)no_omega_0_denominator};
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;
}

void *bath_lorentzian_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
                                                        const grid_volume &gv) const {
  // Compute the number of realnum values required for the polarization data + the bath oscillator
  // data. Each polarization direction is assigned to num_bath oscillators at each spatial grid
  // point
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) num += 2 * gv.ntot() * (1 + num_bath);
  }

  // Allocate memory using the original lorentzian_data structure
  size_t sz = sizeof(lorentzian_data) + sizeof(realnum) * (num - 1);
  lorentzian_data *d = (lorentzian_data *)malloc(sz);
  if (d == NULL) meep::abort("%s:%i:out of memory(%lu)", __FILE__, __LINE__, sz);
  d->sz_data = sz;
  return (void *)d;
}

void bath_lorentzian_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
                                                        realnum dt, const grid_volume &gv,
                                                        void *data) const {
  (void)dt; // unused
  lorentzian_data *d = (lorentzian_data *)data;
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  size_t ntot = d->ntot = gv.ntot();
  realnum *P = d->data;
  realnum *P_prev = d->data + ntot;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) {
      d->P[c][cmp] = P;
      d->P_prev[c][cmp] = P_prev;
      // then the rest space is reserved for the bath oscillators, the pointers for these bath
      // oscillators will be initialized when updating the P field during the dynamics motion
      P += 2 * ntot * (1 + num_bath);
      P_prev += 2 * ntot * (1 + num_bath);
    }
  }

  /*
  master_printf("Using Bath-Lorentzian: num_bath = %d \n", num_bath);
  for (int i = 0; i < num_bath; i++)
  {
    printf("bath freq = %.5E\n", bath_frequencies[i]);
    printf("bath coup = %.5E\n", bath_couplings[i]);
    printf("bath gamma = %.5E\n", bath_gammas[i]);
    printf("bath anharmonicities = %.5E\n", bath_anharmonicities[i]);
  }
  */
  //printf("ntot = %d\n", ntot);
  //printf("size_data = %d\n", sz_data);
  //size_t sz_bath = sizeof(realnum) * 2 * gv.ntot() * num_bath;
  //printf("size_bath = %d\n", sz_bath);

  //master_printf("Bath-Lorentzian param freq = %.5E\n", this->lorentzian_susceptibility::omega_0);
  //master_printf("Bath-Lorentzian param gamma = %.5E\n", this->lorentzian_susceptibility::gamma);
  //printf("conventional Lorentzian param no_omega_0_denominator = %d\n", this->lorentzian_susceptibility::no_omega_0_denominator);

}

void *bath_lorentzian_susceptibility::copy_internal_data(void *data) const {
  lorentzian_data *d = (lorentzian_data *)data;
  if (!d) return 0;
  lorentzian_data *dnew = (lorentzian_data *)malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);
  size_t ntot = d->ntot;
  realnum *P = dnew->data;
  realnum *P_prev = dnew->data + ntot;
  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      dnew->P[c][cmp] = P;
      dnew->P_prev[c][cmp] = P_prev;
      // then the rest space is reserved for the bath oscillators, the pointers for these bath
      // oscillators will be initialized when updating the P field during the dynamics motion
      P += 2 * ntot * (1 + num_bath);
      P_prev += 2 * ntot * (1 + num_bath);
    }
  }
  return (void *)dnew;
}

void bath_lorentzian_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], realnum dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *)P_internal_data;
  const realnum omega2pi = 2 * pi * omega_0, g2pi = gamma * 2 * pi;
  const realnum omega0dtsqr = omega2pi * omega2pi * dt * dt;
  const realnum gamma1inv = 1 / (1 + g2pi * dt / 2), gamma1 = (1 - g2pi * dt / 2);
  const realnum omega0dtsqr_denom = no_omega_0_denominator ? 0 : omega0dtsqr;
  // copy the noise amp from the Noisy Lorentzian susceptibility
  const realnum amp = omega2pi * noise_amp * sqrt(g2pi) * dt * dt / (1 + g2pi * dt / 2); 
  (void)W_prev; // unused;

  // let's define some prefactors necessary for bath_lorentzian calculations
  // first multiply the bath_frequencies and bath_gammas with a factor of 2pi for consistency

  // This piece of code has no problem. It is just very slow. I need to figure out the time consuming steps
  // and then optimize it.
  //typedef std::chrono::high_resolution_clock Clock;
  //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  realnum dtover2 = dt / 2.0;
  realnum twopi = 2 * pi;
  realnum dt2 = dt * dt;

  realnum bathfreq2pi[num_bath], bathgamma2pi[num_bath], bath_couplings2pi[num_bath];
  bool has_anharmonicity = false;
  realnum coeff_a[num_bath], coeff_bplusone[num_bath], coeff_c[num_bath], coeff_d[num_bath], coeff_e[num_bath];
  realnum coeff_ak[num_bath], coeff_bk[num_bath], coeff_dk[num_bath], coeff_ek[num_bath];
  realnum ap = 1.0 + g2pi * dt / 2, prefactor_pnminus = 1.0 - g2pi * dt / 2;

  for (int i = 0; i < num_bath; i++)
  {
    bathfreq2pi[i] = bath_frequencies[i] * twopi;
    bathgamma2pi[i] = bath_gammas[i] * twopi;
    bath_couplings2pi[i] = bath_couplings[i] * twopi;
    // avoid adding the 2pi factor for the bath anharmonicity
    if (abs(bath_anharmonicities[i]) > 1e-20) has_anharmonicity = true; // check if any anharmonicity is present
  
    realnum denom = 1.0 + bathgamma2pi[i] * dtover2;
    realnum ai = (2.0 - bathfreq2pi[i] * bathfreq2pi[i] * dt2) / denom;
    realnum bi = -2.0 / denom;
    realnum ci = bath_couplings2pi[i] * dtover2 / denom;
    realnum di = 1.5 * bath_anharmonicities[i] * dt2 * bathfreq2pi[i] * bathfreq2pi[i] / denom;
    realnum ei = -7.0 / 6.0 * bath_anharmonicities[i] * bath_anharmonicities[i] * dt2 * bathfreq2pi[i] * bathfreq2pi[i] / denom;

    coeff_a[i] = ai;
    coeff_c[i] = ci;
    coeff_d[i] = di;
    coeff_e[i] = ei;
    coeff_ak[i] = dtover2 * ai * bath_couplings2pi[i];
    coeff_bk[i] = dtover2 * bi * bath_couplings2pi[i];
    coeff_dk[i] = dtover2 * di * bath_couplings2pi[i];
    coeff_ek[i] = dtover2 * ei * bath_couplings2pi[i];
    coeff_bplusone[i] = bi + 1.0;
    
    realnum coupling_term = dtover2 * bath_couplings2pi[i] * ci;
    ap += coupling_term;
    prefactor_pnminus -= coupling_term;
  }
  realnum apinv = 1.0 / ap;


  //std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  //std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t2 - t1);
  //std::cout << "Update Bath Lorentz medium for one time step before for loop takes " << time_span.count() << " seconds." << std::endl;

  // TODO: add back lorentzian_unstable(omega_0, gamma, dt) if we can improve the stability test
  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
      if (w && s) {
        realnum *p = d->P[c][cmp], *pp = d->P_prev[c][cmp];
        // also create pointers for the bath oscillators

        /*
        realnum *p_bath[num_bath];
        realnum *pp_bath[num_bath];
        for (int k = 0; k < num_bath; k++)
        {
            p_bath[k] = pp + d->ntot + d->ntot * k * 2;
            pp_bath[k] = pp + d->ntot + d->ntot * (k * 2 + 1);
        }
        */
        
        // new code for improving the performance when iterating the bath field
        size_t ntot = d->ntot;
        realnum *bath_data_start = pp + ntot;
        realnum *p_bath_base = bath_data_start;
        realnum *pp_bath_base = bath_data_start + (num_bath * ntot);

        // directions/strides for offdiagonal terms, similar to update_eh
        const direction d = component_direction(c);
        const ptrdiff_t is = gv.stride(d) * (is_magnetic(c) ? -1 : +1);
        direction d1 = cycle_direction(gv.dim, d, 1);
        component c1 = direction_component(c, d1);
        ptrdiff_t is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
        const realnum *w1 = W[c1][cmp];
        const realnum *s1 = w1 ? sigma[c][d1] : NULL;
        direction d2 = cycle_direction(gv.dim, d, 2);
        component c2 = direction_component(c, d2);
        ptrdiff_t is2 = gv.stride(d2) * (is_magnetic(c) ? -1 : +1);
        const realnum *w2 = W[c2][cmp];
        const realnum *s2 = w2 ? sigma[c][d2] : NULL;
        
        // for 3x3 or 2x2 anisotropic systems, we keep the original code
        if (s2 && !s1) { // make s1 the non-NULL one if possible
          SWAP(direction, d1, d2);
          SWAP(component, c1, c2);
          SWAP(ptrdiff_t, is1, is2);
          SWAP(const realnum *, w1, w2);
          SWAP(const realnum *, s1, s2);
        }
        if (s1 && s2) { // 3x3 anisotropic
          PLOOP_OVER_VOL_OWNED(gv, c, i) {
            // s[i] != 0 check is a bit of a hack to work around
            // some instabilities that occur near the boundaries
            // of materials; see PR #666
            if (s[i] != 0) {
              realnum pcur = p[i];
              p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] +
                                  omega0dtsqr * (s[i] * w[i] + OFFDIAG(s1, w1, is1, is) +
                                                 OFFDIAG(s2, w2, is2, is)));
              pp[i] = pcur;
            }
          }
        }
        else if (s1) { // 2x2 anisotropic
          PLOOP_OVER_VOL_OWNED(gv, c, i) {
            if (s[i] != 0) { // see above
              realnum pcur = p[i];
              p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] +
                                  omega0dtsqr * (s[i] * w[i] + OFFDIAG(s1, w1, is1, is)));
              pp[i] = pcur;
            }
          }
        }
        else { // isotropic
          // We only implement the bath-Lorentz model for isotropic systems
          // TODO
          PLOOP_OVER_VOL_OWNED(gv, c, i) {
            // This is the equation of motion for an independent Lorentz oscillator, the damping term comes from the gamma1
            realnum pcur = p[i];
            // new code for improving the performance when iterating the bath field
            realnum *p_bath_i = p_bath_base + i * num_bath;
            realnum *pp_bath_i = pp_bath_base + i * num_bath;
            //std::vector<realnum> pbathcur, pbathpre;
            //for(size_t k = 0; k< num_bath; k++) 
            //{
            //  pbathcur.push_back(p_bath[k][i]);
            //  pbathpre.push_back(pp_bath[k][i]);
            //}
            realnum pbathcur[num_bath];
            realnum pbathcur2[num_bath];
            realnum pbathcur3[num_bath];
            realnum sum_kiaiYi_cur = 0;
            realnum sum_kibiYi_pre = 0;
            realnum sum_kidiYi3_cur = 0;
            // #pragma omp simd reduction(+:sum_kiaiYi_cur, sum_kibiYi_pre, sum_kidiYi3_cur)
            for(int k = 0; k< num_bath; k++) 
            {
              /*
              pbathcur[k] = p_bath[k][i];
              sum_kiaiYi_cur += coeff_ak[k] * pbathcur[k];
              sum_kibiYi_pre += coeff_bk[k] * pp_bath[k][i];
              if (has_anharmonicity)
                sum_kidiYi3_cur += coeff_dk[k] * pbathcur[k] * pbathcur[k] + coeff_ek[k] * pbathcur[k] * pbathcur[k] * pbathcur[k];
              */

              // new code for improving the performance when iterating the bath field
              realnum pbc_k = p_bath_i[k];
              pbathcur[k] = pbc_k; 
              pbathcur2[k] = pbc_k * pbc_k;
              pbathcur3[k] = pbathcur2[k] * pbc_k;
              sum_kiaiYi_cur += coeff_ak[k] * pbc_k;
              sum_kibiYi_pre += coeff_bk[k] * pp_bath_i[k];
              if (has_anharmonicity)
                sum_kidiYi3_cur += coeff_dk[k] * pbathcur2[k] + coeff_ek[k] * pbathcur3[k];
            }

            // precompute some important quantities
            //realnum sum_kiaiYi_cur = dt / 2.0 * std::inner_product(coeff_ak.begin(), coeff_ak.end(), pbathcur.begin(), 0.0);
            //realnum sum_kibiYi_pre = dt / 2.0 * std::inner_product(coeff_bk.begin(), coeff_bk.end(), pbathpre.begin(), 0.0);
            // update P to the next time step
            p[i] = apinv * (pcur * (2 - omega0dtsqr_denom) - prefactor_pnminus * pp[i] + omega0dtsqr * (s[i] * w[i]) - sum_kiaiYi_cur - sum_kibiYi_pre - sum_kidiYi3_cur);
            // update bath coordinates to the next time step
            realnum p_pp_diff = p[i] - pp[i];
            double gaussian_random_amp = amp * sqrt(s[i]);
            /*
            for (int k = 0; k < num_bath; k++)
            {
              //p_bath[k][i] = coeff_a[k] * pbathcur[k] + (coeff_b[k] + 1.0) * pbathpre[k] + coeff_c[k] * (p[i] - pp[i]);
              double anharmonicity_term = 0.0;
              if (has_anharmonicity)
                 anharmonicity_term = coeff_d[k] * pbathcur[k] * pbathcur[k] + coeff_e[k] * pbathcur[k] * pbathcur[k] * pbathcur[k];
              p_bath[k][i] = coeff_a[k] * pbathcur[k] + coeff_bplusone[k] * pp_bath[k][i] + coeff_c[k] * p_pp_diff + anharmonicity_term;
              // consider to add a noisy term to account for the thermal fluctuations of the bath oscillators
              if (noise_amp > 1e-10)
                p_bath[k][i] += gaussian_random(0, gaussian_random_amp);
              // reset the previous values
              pp_bath[k][i] = pbathcur[k];
            }
            */

            // new code for improving the performance when iterating the bath field
            // #pragma omp simd
            // reduce the cost for Gaussian random number generation by initializing the seed only once
            // std::mt19937_64 random;
            // replace mt19937_64 with a faster pcg seed
            //pcg_extras::seed_seq_from<std::random_device> seed_source;
            //pcg32_fast random(seed_source);
            //cxx::ziggurat_normal_distribution<double> normal_distr(0, gaussian_random_amp);

            if (!has_anharmonicity)
            {
              if (noise_amp <= 1e-10)
              {
                // no anharmonicity and no noise, nothing to do
                for (int k = 0; k < num_bath; k++)
            {
              //p_bath[k][i] = coeff_a[k] * pbathcur[k] + (coeff_b[k] + 1.0) * pbathpre[k] + coeff_c[k] * (p[i] - pp[i]);
              //double anharmonicity_term = 0.0;
              //if (has_anharmonicity)
              //  anharmonicity_term = coeff_d[k] * pbathcur2[k] + coeff_e[k] * pbathcur3[k];
              p_bath_i[k] = coeff_a[k] * pbathcur[k] + coeff_bplusone[k] * pp_bath_i[k] + coeff_c[k] * p_pp_diff; // + anharmonicity_term;
              // consider to add a noisy term to account for the thermal fluctuations of the bath oscillators
              //if (noise_amp > 1e-10)
                //p_bath_i[k] += gaussian_random(0, gaussian_random_amp);
                // the zigguart gaussian number generator is x5 faster than the above line 
              //  p_bath_i[k] += normal_distr(random);
              // reset the previous values
              pp_bath_i[k] = pbathcur[k];
            }
            
              }
              else
              {
                pcg_extras::seed_seq_from<std::random_device> seed_source;
                pcg32_fast random(seed_source);
                cxx::ziggurat_normal_distribution<double> normal_distr(0, gaussian_random_amp);
                for (int k = 0; k < num_bath; k++)
            {
              //p_bath[k][i] = coeff_a[k] * pbathcur[k] + (coeff_b[k] + 1.0) * pbathpre[k] + coeff_c[k] * (p[i] - pp[i]);
              //double anharmonicity_term = 0.0;
              //if (has_anharmonicity)
              //  anharmonicity_term = coeff_d[k] * pbathcur2[k] + coeff_e[k] * pbathcur3[k];
              p_bath_i[k] = coeff_a[k] * pbathcur[k] + coeff_bplusone[k] * pp_bath_i[k] + coeff_c[k] * p_pp_diff; // + anharmonicity_term;
              // consider to add a noisy term to account for the thermal fluctuations of the bath oscillators
              if (noise_amp > 1e-10)
                //p_bath_i[k] += gaussian_random(0, gaussian_random_amp);
                // the zigguart gaussian number generator is x5 faster than the above line 
                p_bath_i[k] += normal_distr(random);
              // reset the previous values
              pp_bath_i[k] = pbathcur[k];
            }
              }
            }
            else
            {
              if (noise_amp <= 1e-10)
              {
                for (int k = 0; k < num_bath; k++)
            {
              //p_bath[k][i] = coeff_a[k] * pbathcur[k] + (coeff_b[k] + 1.0) * pbathpre[k] + coeff_c[k] * (p[i] - pp[i]);
              double anharmonicity_term = 0.0;
              if (has_anharmonicity)
                anharmonicity_term = coeff_d[k] * pbathcur2[k] + coeff_e[k] * pbathcur3[k];
              p_bath_i[k] = coeff_a[k] * pbathcur[k] + coeff_bplusone[k] * pp_bath_i[k] + coeff_c[k] * p_pp_diff + anharmonicity_term;
              // consider to add a noisy term to account for the thermal fluctuations of the bath oscillators
              //if (noise_amp > 1e-10)
                //p_bath_i[k] += gaussian_random(0, gaussian_random_amp);
                // the zigguart gaussian number generator is x5 faster than the above line 
              //  p_bath_i[k] += normal_distr(random);
              // reset the previous values
              pp_bath_i[k] = pbathcur[k];
            }
              }
              else
              {
                pcg_extras::seed_seq_from<std::random_device> seed_source;
                pcg32_fast random(seed_source);
                cxx::ziggurat_normal_distribution<double> normal_distr(0, gaussian_random_amp);
                for (int k = 0; k < num_bath; k++)
            {
              //p_bath[k][i] = coeff_a[k] * pbathcur[k] + (coeff_b[k] + 1.0) * pbathpre[k] + coeff_c[k] * (p[i] - pp[i]);
              double anharmonicity_term = 0.0;
              if (has_anharmonicity)
                anharmonicity_term = coeff_d[k] * pbathcur2[k] + coeff_e[k] * pbathcur3[k];
              p_bath_i[k] = coeff_a[k] * pbathcur[k] + coeff_bplusone[k] * pp_bath_i[k] + coeff_c[k] * p_pp_diff + anharmonicity_term;
              // consider to add a noisy term to account for the thermal fluctuations of the bath oscillators
              if (noise_amp > 1e-10)
                //p_bath_i[k] += gaussian_random(0, gaussian_random_amp);
                // the zigguart gaussian number generator is x5 faster than the above line 
                p_bath_i[k] += normal_distr(random);
              // reset the previous values
              pp_bath_i[k] = pbathcur[k];
            }
              }
            }

            pp[i] = pcur;

          }
        }
      }
    }
  }

  //std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
  //std::chrono::duration<double> time_span3 = std::chrono::duration_cast< std::chrono::duration<double> >(t3 - t2);
  //std::cout << "Update Bath Lorentz medium for loop takes " << time_span3.count() << " seconds." << std::endl;
}

void bath_lorentzian_susceptibility::dump_params(h5file *h5f, size_t *start) {
  // Total parameters: 5 base + 1 for num_bath + 4 per bath oscillator.
  size_t num_params = 7 + num_bath * 4;
  size_t params_dims[1] = {num_params};

  // Allocate a dynamic array to hold all parameters.
  realnum *params_data = new realnum[num_params];

  // Fill in the base parameters.
  params_data[0] = num_params - 1;
  params_data[1] = (realnum)get_id();
  params_data[2] = omega_0;
  params_data[3] = gamma;
  params_data[4] = (realnum)no_omega_0_denominator;

  // Add the number of bath oscillators.
  params_data[5] = (realnum)num_bath;
  params_data[6] = (realnum)noise_amp;

  // Fill in the bath oscillator parameters.
  size_t index = 7;
  for (int i = 0; i < num_bath; ++i) {
    params_data[index++] = bath_frequencies[i];
    params_data[index++] = bath_couplings[i];
    params_data[index++] = bath_gammas[i];
    params_data[index++] = bath_anharmonicities[i];
  }

  // Write the chunk.
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;

  // Clean up the dynamic array.
  delete[] params_data;
}

mxl_socket_susceptibility::mxl_socket_susceptibility(realnum rescaling_factor,
                                                     realnum time_units_fs, realnum timeout,
                                                     const char *host, int port,
                                                     const char *label,
                                                     bool real_field_only)
    : rescaling_factor(rescaling_factor), time_units_fs(time_units_fs), timeout(timeout),
      host(host ? host : "127.0.0.1"), port(port), label(label ? label : ""),
      real_field_only(real_field_only) {
  mxl_socket_susceptibility_present = true;
  mxl_assert_little_endian();
  if (rescaling_factor < 0.0)
    meep::abort("MXLSocketSusceptibility rescaling_factor must be nonnegative.");
  if (time_units_fs <= 0.0)
    meep::abort("MXLSocketSusceptibility time_units_fs must be positive.");
  if (timeout < 0.0) meep::abort("MXLSocketSusceptibility timeout must be nonnegative.");
  if (this->host.empty()) meep::abort("MXLSocketSusceptibility host must be nonempty.");
  if (port <= 0 || port > 65535)
    meep::abort("MXLSocketSusceptibility port must be in the range [1, 65535].");
  if (this->label.find('\n') != std::string::npos ||
      this->label.find('\r') != std::string::npos)
    meep::abort("MXLSocketSusceptibility label must not contain a newline.");
}

bool mxl_socket_susceptibility::needs_P(component c, int cmp,
                                        realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  if (!is_electric(c)) return false;
  direction d = component_direction(c);
  /* In complex-field Meep runs, boundary exchange expects valid real and
     imaginary polarization storage.  In real_field_only mode P[c][1] remains
     allocated but is left at zero; otherwise independent socket molecules
     drive and update it from W[c][1]. */
  return !trivial_sigma[c][d] && W[c][cmp];
}

bool mxl_socket_susceptibility::needs_W_notowned(component c,
                                                 realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  (void)c;
  (void)W;
  return false;
}

void *mxl_socket_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
                                                   const grid_volume &gv) const {
  (void)W;
  (void)gv;
  return new mxl_socket_data();
}

void mxl_socket_susceptibility::delete_internal_data(void *data) const {
  delete (mxl_socket_data *)data;
}

void mxl_socket_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
                                                   realnum dt, const grid_volume &gv,
                                                   void *data) const {
  (void)dt;
  mxl_socket_data *d = (mxl_socket_data *)data;
  d->ntot = gv.ntot();

  FOR_COMPONENTS(c) DOCMP2 {
    d->P_store[c][cmp].clear();
    d->P[c][cmp] = NULL;
    if (needs_P(c, cmp, W)) d->P_store[c][cmp].assign(d->ntot, 0.0);
  }
  d->reset_ptrs();

  std::vector<char> active(d->ntot, 0);
  /* The sigma arrays provide the material active-site mask here, not the
     oscillator-strength semantics used by Lorentzian media. */
  FOR_COMPONENTS(c) {
    if (!is_electric(c) || !gv.has_field(c)) continue;
    direction dc = component_direction(c);
    const realnum *s = sigma[c][dc];
    if (!s || (!W[c][0] && !W[c][1])) continue;
    PLOOP_OVER_VOL_OWNED(gv, c, i) {
      if (s[i] != 0.0) {
        int axis = mxl_amp_axis(c); 
        // 1 << axis (axis = 0,1,2 for xyz) is a bitmask to obtain 001, 010, 100 (1, 2, 4)
        int component_bit = 1 << axis; 
        active[(size_t)i] = active[(size_t)i] | component_bit;
      }

    }
  }

  d->active_indices.clear();
  d->comp_mask.clear();
  for (size_t i = 0; i < active.size(); ++i)
    if (active[i]) {
      d->active_indices.push_back(i);
      d->comp_mask.push_back(active[i]);
    }

  d->drive_cmps.clear();
  d->drive_cmps.push_back(0);
  if (!real_field_only) {
    bool has_imag_fields = false;
    FOR_COMPONENTS(c) {
      if (is_electric(c) && gv.has_field(c) && W[c][1]) {
        has_imag_fields = true;
        break;
      }
    }
    if (has_imag_fields) d->drive_cmps.push_back(1);
  }

  const size_t required_driver_count = d->active_indices.size() * d->drive_cmps.size();
  if (required_driver_count >= (size_t)mxl_molecules_per_chunk_id_block)
    meep::abort("MXLSocketSusceptibility required socket driver count %zu exceeds per-chunk id "
                "block %d.",
                required_driver_count, mxl_molecules_per_chunk_id_block);

  /* Rank/chunk id blocks keep molecule ids deterministic for a fixed chunk
     decomposition. */
  int chunk_ordinal = mxl_next_chunk_ordinal++;
  if (chunk_ordinal >= mxl_chunks_per_rank_id_block)
    meep::abort("MXLSocketSusceptibility rank %d has too many internal chunks (%d >= %d).",
                my_rank(), chunk_ordinal + 1, mxl_chunks_per_rank_id_block);

  const long long rank_chunk =
      ((long long)my_rank()) * mxl_chunks_per_rank_id_block + chunk_ordinal;
  const long long base_id_long = rank_chunk * mxl_molecules_per_chunk_id_block;
  if (base_id_long + (long long)required_driver_count > mxl_max_molecule_id)
    meep::abort("MXLSocketSusceptibility molecule ids exceed the int32 protocol range.");

  int base_id = (int)base_id_long;
  d->molecule_ids.resize(required_driver_count);
  for (size_t i = 0; i < d->molecule_ids.size(); ++i)
    d->molecule_ids[i] = base_id + (int)i;

  mxl_local_mapping_records.reserve(mxl_local_mapping_records.size() + required_driver_count);
  for (size_t icmp = 0; icmp < d->drive_cmps.size(); ++icmp) {
    int cmp = d->drive_cmps[icmp];
    for (size_t isite = 0; isite < d->active_indices.size(); ++isite) {
      size_t imol = icmp * d->active_indices.size() + isite;
      size_t idx = d->active_indices[isite];
      vec r = gv.loc(Centered, (ptrdiff_t)idx);
      mxl_socket_mapping_record rec;
      rec.molecule_id = d->molecule_ids[imol];
      rec.chunk_ordinal = chunk_ordinal;
      rec.drive_cmp = cmp;
      rec.site_ordinal = isite;
      rec.local_index = idx;
      rec.x_or_r = mxl_x_or_r(r);
      rec.y = mxl_y(r);
      rec.z = mxl_z(r);
      rec.label = label;
      rec.host = host;
      rec.port = port;
      mxl_local_mapping_records.push_back(rec);
    }
  }

  mxl_local_active_site_count += d->active_indices.size();
  mxl_local_required_driver_count += required_driver_count;
  if (!d->active_indices.empty() && d->drive_cmps.size() > 1)
    mxl_socket_imag_field_coupling_active = true;

  d->efields_au.assign(3 * required_driver_count, 0.0);
  d->amps_au.assign(3 * required_driver_count, 0.0);
  d->initialized = false;
}

void *mxl_socket_susceptibility::copy_internal_data(void *data) const {
  mxl_socket_data *d = (mxl_socket_data *)data;
  if (!d) return NULL;
  return new mxl_socket_data(*d);
}

void mxl_socket_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], realnum dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  (void)W_prev;
  mxl_socket_data *d = (mxl_socket_data *)P_internal_data;
  if (!d || d->active_indices.empty() || rescaling_factor == 0.0) return;

  const double efield_factor =
      mxl_efield_mu_to_au_prefactor * rescaling_factor / (time_units_fs * time_units_fs);
  const size_t nsites = d->active_indices.size();

  component comps[3];
  mxl_field_components(gv.dim, comps);
  for (size_t icmp = 0; icmp < d->drive_cmps.size(); ++icmp) {
    int cmp = d->drive_cmps[icmp];
    for (size_t isite = 0; isite < nsites; ++isite) {
      size_t imol = icmp * nsites + isite;
      size_t idx = d->active_indices[isite];
      for (int axis = 0; axis < 3; ++axis) {
        component c = comps[axis];
        d->efields_au[3 * imol + axis] =
            ((d->comp_mask[isite] & (1 << axis)) && W[c][cmp])
                ? efield_factor * W[c][cmp][idx]
                : 0.0;
      }
    }
  }

  d->client.ensure_connected(host, port, timeout);
  if (!d->initialized) {
    double dt_au = dt * time_units_fs * mxl_fs_to_au;
    d->client.send_init(d->molecule_ids, dt_au, rescaling_factor, time_units_fs, timeout,
                        my_rank(), mxl_expected_total_required_driver_count);
    d->initialized = true;
  }
  d->client.step(d->molecule_ids, d->efields_au, d->amps_au);

  const int dim_power = mxl_dim_power(gv.dim);
  double cell_measure = 1.0;
  for (int i = 0; i < dim_power; ++i)
    cell_measure *= gv.inva;
  /* Convert molecular dmu/dt into a polarization-density increment. */
  const double pdot_scale = mxl_source_amp_au_to_mu * rescaling_factor / cell_measure;

  FOR_COMPONENTS(c) {
    if (!is_electric(c) || !gv.has_field(c)) continue;
    direction dc = component_direction(c);
    const realnum *s = sigma[c][dc];
    if (!s) continue;
    int axis = mxl_amp_axis(c);
    for (size_t icmp = 0; icmp < d->drive_cmps.size(); ++icmp) {
      int cmp = d->drive_cmps[icmp];
      realnum *p = d->P[c][cmp];
      if (!p) continue;
      for (size_t isite = 0; isite < nsites; ++isite) {
        size_t imol = icmp * nsites + isite;
        size_t idx = d->active_indices[isite];
        if (d->comp_mask[isite] & (1 << axis))
          p[idx] += (realnum)(dt * pdot_scale * s[idx] * d->amps_au[3 * imol + axis]);
      }
    }
  }
}

void mxl_socket_susceptibility::subtract_P(field_type ft,
                                           realnum *f_minus_p[NUM_FIELD_COMPONENTS][2],
                                           void *P_internal_data) const {
  mxl_socket_data *d = (mxl_socket_data *)P_internal_data;
  if (!d) return;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff;
  FOR_FT_COMPONENTS(ft, ec) DOCMP2 {
    if (d->P[ec][cmp]) {
      component dc = field_type_component(ft2, ec);
      if (f_minus_p[dc][cmp]) {
        realnum *p = d->P[ec][cmp];
        realnum *fmp = f_minus_p[dc][cmp];
        for (size_t i = 0; i < d->ntot; ++i)
          fmp[i] -= p[i];
      }
    }
  }
}

int mxl_socket_susceptibility::num_cinternal_notowned_needed(component c,
                                                             void *P_internal_data) const {
  mxl_socket_data *d = (mxl_socket_data *)P_internal_data;
  return d && d->P[c][0] ? 1 : 0;
}

realnum *mxl_socket_susceptibility::cinternal_notowned_ptr(int inotowned, component c, int cmp,
                                                           int n, void *P_internal_data) const {
  (void)inotowned;
  mxl_socket_data *d = (mxl_socket_data *)P_internal_data;
  if (!d || !d->P[c][cmp]) return NULL;
  return d->P[c][cmp] + n;
}

void mxl_socket_susceptibility::dump_params(h5file *h5f, size_t *start) {
  (void)h5f;
  (void)start;
  meep::abort("HDF5 dumping of MXLSocketSusceptibility is not supported.");
}

int mxl_socket_susceptibility::get_num_params() {
  meep::abort("HDF5 dumping of MXLSocketSusceptibility is not supported.");
  return 0;
}

gyrotropic_susceptibility::gyrotropic_susceptibility(const vec &bias, realnum omega_0,
                                                     realnum gamma, realnum alpha,
                                                     gyrotropy_model model)
    : omega_0(omega_0), gamma(gamma), alpha(alpha), model(model) {
  // Precalculate g_{ij} = sum_k epsilon_{ijk} b_k, used in update_P.
  // Ignore |b| for Landau-Lifshitz-Gilbert gyrotropy model.
  const vec b = (model == GYROTROPIC_SATURATED) ? bias / abs(bias) : bias;
  memset(gyro_tensor, 0, 9 * sizeof(realnum));
  gyro_tensor[X][Y] = b.z();
  gyro_tensor[Y][X] = -b.z();
  gyro_tensor[Y][Z] = b.x();
  gyro_tensor[Z][Y] = -b.x();
  gyro_tensor[Z][X] = b.y();
  gyro_tensor[X][Z] = -b.y();
}

/* To implement gyrotropic susceptibilities, we track three
   polarization components (e.g. Px, Py, Pz) on EACH of the Yee cell's
   three driving field positions (e.g., Ex, Ey, and Ez), i.e. 9
   numbers per cell.  This takes 3X the memory and runtime compared to
   Lorentzian susceptibility.  The advantage is that during update_P,
   we can directly access the value of P at each update point without
   averaging.  */

typedef struct {
  size_t sz_data;
  size_t ntot;
  realnum *P[NUM_FIELD_COMPONENTS][2][3];
  realnum *P_prev[NUM_FIELD_COMPONENTS][2][3];
  realnum data[1];
} gyrotropy_data;

void *gyrotropic_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
                                                   const grid_volume &gv) const {
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) num += 6 * gv.ntot();
  }
  size_t sz = sizeof(gyrotropy_data) + sizeof(realnum) * (num - 1);
  gyrotropy_data *d = (gyrotropy_data *)malloc(sz);
  if (d == NULL) meep::abort("%s:%i:out of memory(%lu)", __FILE__, __LINE__, sz);
  d->sz_data = sz;
  return (void *)d;
}

void gyrotropic_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], realnum dt,
                                                   const grid_volume &gv, void *data) const {
  (void)dt; // unused
  gyrotropy_data *d = (gyrotropy_data *)data;
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  d->ntot = gv.ntot();
  realnum *p = d->data;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) {
      for (int dd = X; dd < R; dd++) {
        d->P[c][cmp][dd] = p;
        p += d->ntot;
        d->P_prev[c][cmp][dd] = p;
        p += d->ntot;
      }
    }
  }
}

void *gyrotropic_susceptibility::copy_internal_data(void *data) const {
  gyrotropy_data *d = (gyrotropy_data *)data;
  if (!d) return 0;
  gyrotropy_data *dnew = (gyrotropy_data *)malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);
  realnum *p = dnew->data;
  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp][0]) {
      for (int dd = X; dd < R; dd++) {
        dnew->P[c][cmp][dd] = p;
        p += d->ntot;
        dnew->P_prev[c][cmp][dd] = p;
        p += d->ntot;
      }
    }
  }
  return (void *)dnew;
}

bool gyrotropic_susceptibility::needs_P(component c, int cmp,
                                        realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  if (!is_electric(c) && !is_magnetic(c)) return false;
  direction d0 = component_direction(c);
  return (d0 == X || d0 == Y || d0 == Z) && sigma[c][d0] && W[c][cmp];
}

// Similar to the OFFDIAG macro, but without averaging sigma.
#define OFFDIAGW(g, sx, s) (0.25 * (g[i] + g[i - sx] + g[i + s] + g[i + s - sx]))

void gyrotropic_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], realnum dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  gyrotropy_data *d = (gyrotropy_data *)P_internal_data;
  const realnum omega2pidt = 2 * pi * omega_0 * dt;
  const realnum g2pidt = 2 * pi * gamma * dt;
  (void)W_prev; // unused;

  switch (model) {
    case GYROTROPIC_LORENTZIAN:
    case GYROTROPIC_DRUDE: {
      const realnum omega0dtsqr = omega2pidt * omega2pidt;
      const realnum gamma1 = (1 - g2pidt / 2);
      const realnum diag = 2 - (model == GYROTROPIC_DRUDE ? 0 : omega0dtsqr);
      const realnum pt = pi * dt;

      // Precalculate 3x3 matrix inverse, exploiting skew symmetry
      const realnum gd = (1 + g2pidt / 2);
      const realnum gx = pt * gyro_tensor[Y][Z];
      const realnum gy = pt * gyro_tensor[Z][X];
      const realnum gz = pt * gyro_tensor[X][Y];
      const realnum invdet = 1.0 / gd / (gd * gd + gx * gx + gy * gy + gz * gz);
      const realnum inv[3][3] = {{invdet * (gd * gd + gx * gx), invdet * (gx * gy + gd * gz),
                                  invdet * (gx * gz - gd * gy)},
                                 {invdet * (gy * gx - gd * gz), invdet * (gd * gd + gy * gy),
                                  invdet * (gy * gz + gd * gx)},
                                 {invdet * (gz * gx + gd * gy), invdet * (gz * gy - gd * gx),
                                  invdet * (gd * gd + gz * gz)}};

      FOR_COMPONENTS(c) DOCMP2 {
        if (d->P[c][cmp][0]) {
          const direction d0 = component_direction(c);
          const realnum *w0 = W[c][cmp], *s = sigma[c][d0];

          if (!w0 || !s || (d0 != X && d0 != Y && d0 != Z))
            meep::abort("gyrotropic media require 3D Cartesian fields\n");

          const direction d1 = cycle_direction(gv.dim, d0, 1);
          const direction d2 = cycle_direction(gv.dim, d0, 2);
          const realnum *w1 = W[direction_component(c, d1)][cmp];
          const realnum *w2 = W[direction_component(c, d2)][cmp];
          realnum *p0 = d->P[c][cmp][d0], *pp0 = d->P_prev[c][cmp][d0];
          realnum *p1 = d->P[c][cmp][d1], *pp1 = d->P_prev[c][cmp][d1];
          realnum *p2 = d->P[c][cmp][d2], *pp2 = d->P_prev[c][cmp][d2];
          const ptrdiff_t is = gv.stride(d0) * (is_magnetic(c) ? -1 : +1);
          const ptrdiff_t is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
          const ptrdiff_t is2 = gv.stride(d2) * (is_magnetic(c) ? -1 : +1);
          realnum r0, r1, r2;

          if (!pp1 || !pp2) meep::abort("gyrotropic media require 3D Cartesian fields\n");
          if (sigma[c][d1] || sigma[c][d2])
            meep::abort("gyrotropic media do not support anisotropic sigma\n");

          LOOP_OVER_VOL_OWNED(gv, c, i) {
            r0 = diag * p0[i] - gamma1 * pp0[i] + omega0dtsqr * s[i] * w0[i] -
                 pt * gyro_tensor[d0][d1] * pp1[i] - pt * gyro_tensor[d0][d2] * pp2[i];
            r1 = diag * p1[i] - gamma1 * pp1[i] +
                 (w1 ? omega0dtsqr * s[i] * OFFDIAGW(w1, is1, is) : 0) -
                 pt * gyro_tensor[d1][d0] * pp0[i] - pt * gyro_tensor[d1][d2] * pp2[i];
            r2 = diag * p2[i] - gamma1 * pp2[i] +
                 (w2 ? omega0dtsqr * s[i] * OFFDIAGW(w2, is2, is) : 0) -
                 pt * gyro_tensor[d2][d1] * pp1[i] - pt * gyro_tensor[d2][d0] * pp0[i];

            pp0[i] = p0[i];
            pp1[i] = p1[i];
            pp2[i] = p2[i];
            p0[i] = inv[d0][d0] * r0 + inv[d0][d1] * r1 + inv[d0][d2] * r2;
            p1[i] = inv[d1][d0] * r0 + inv[d1][d1] * r1 + inv[d1][d2] * r2;
            p2[i] = inv[d2][d0] * r0 + inv[d2][d1] * r1 + inv[d2][d2] * r2;
          }
        }
      }
    } break;

    case GYROTROPIC_SATURATED: {
      const realnum dt2pi = 2 * pi * dt;

      // Precalculate 3x3 matrix inverse, exploiting skew symmetry
      const realnum gd = 0.5;
      const realnum gx = -0.5 * alpha * gyro_tensor[Y][Z];
      const realnum gy = -0.5 * alpha * gyro_tensor[Z][X];
      const realnum gz = -0.5 * alpha * gyro_tensor[X][Y];
      const realnum invdet = 1.0 / gd / (gd * gd + gx * gx + gy * gy + gz * gz);
      const realnum inv[3][3] = {{invdet * (gd * gd + gx * gx), invdet * (gx * gy + gd * gz),
                                  invdet * (gx * gz - gd * gy)},
                                 {invdet * (gy * gx - gd * gz), invdet * (gd * gd + gy * gy),
                                  invdet * (gy * gz + gd * gx)},
                                 {invdet * (gz * gx + gd * gy), invdet * (gz * gy - gd * gx),
                                  invdet * (gd * gd + gz * gz)}};

      FOR_COMPONENTS(c) DOCMP2 {
        if (d->P[c][cmp][0]) {
          const direction d0 = component_direction(c);
          const realnum *w0 = W[c][cmp], *s = sigma[c][d0];

          if (!w0 || !s || (d0 != X && d0 != Y && d0 != Z))
            meep::abort("gyrotropic media require 3D Cartesian fields\n");

          const direction d1 = cycle_direction(gv.dim, d0, 1);
          const direction d2 = cycle_direction(gv.dim, d0, 2);
          const realnum *w1 = W[direction_component(c, d1)][cmp];
          const realnum *w2 = W[direction_component(c, d2)][cmp];
          realnum *p0 = d->P[c][cmp][d0], *pp0 = d->P_prev[c][cmp][d0];
          realnum *p1 = d->P[c][cmp][d1], *pp1 = d->P_prev[c][cmp][d1];
          realnum *p2 = d->P[c][cmp][d2], *pp2 = d->P_prev[c][cmp][d2];
          const ptrdiff_t is = gv.stride(d0) * (is_magnetic(c) ? -1 : +1);
          const ptrdiff_t is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
          const ptrdiff_t is2 = gv.stride(d2) * (is_magnetic(c) ? -1 : +1);
          realnum r0, r1, r2, q0, q1, q2;

          if (!pp1 || !pp2) meep::abort("gyrotropic media require 3D Cartesian fields\n");
          if (sigma[c][d1] || sigma[c][d2])
            meep::abort("gyrotropic media do not support anisotropic sigma\n");

          LOOP_OVER_VOL_OWNED(gv, c, i) {
            q0 = -omega2pidt * p0[i] + 0.5 * alpha * pp0[i] + dt2pi * s[i] * w0[i];
            q1 = -omega2pidt * p1[i] + 0.5 * alpha * pp1[i] +
                 dt2pi * s[i] * (w1 ? OFFDIAGW(w1, is1, is) : 0);
            q2 = -omega2pidt * p2[i] + 0.5 * alpha * pp2[i] +
                 dt2pi * s[i] * (w2 ? OFFDIAGW(w2, is2, is) : 0);

            r0 =
                0.5 * pp0[i] - g2pidt * p0[i] + gyro_tensor[d0][d1] * q1 + gyro_tensor[d0][d2] * q2;
            r1 =
                0.5 * pp1[i] - g2pidt * p1[i] + gyro_tensor[d1][d2] * q2 + gyro_tensor[d1][d0] * q0;
            r2 =
                0.5 * pp2[i] - g2pidt * p2[i] + gyro_tensor[d2][d0] * q0 + gyro_tensor[d2][d1] * q1;

            pp0[i] = p0[i];
            pp1[i] = p1[i];
            pp2[i] = p2[i];
            p0[i] = inv[d0][d0] * r0 + inv[d0][d1] * r1 + inv[d0][d2] * r2;
            p1[i] = inv[d1][d0] * r0 + inv[d1][d1] * r1 + inv[d1][d2] * r2;
            p2[i] = inv[d2][d0] * r0 + inv[d2][d1] * r1 + inv[d2][d2] * r2;
          }
        }
      }
    } break;
  }
}

void gyrotropic_susceptibility::subtract_P(field_type ft,
                                           realnum *f_minus_p[NUM_FIELD_COMPONENTS][2],
                                           void *P_internal_data) const {
  gyrotropy_data *d = (gyrotropy_data *)P_internal_data;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  size_t ntot = d->ntot;
  FOR_FT_COMPONENTS(ft, ec) DOCMP2 {
    if (d->P[ec][cmp][0]) {
      component dc = field_type_component(ft2, ec);
      if (f_minus_p[dc][cmp]) {
        realnum *p = d->P[ec][cmp][component_direction(ec)];
        realnum *fmp = f_minus_p[dc][cmp];
        for (size_t i = 0; i < ntot; ++i)
          fmp[i] -= p[i];
      }
    }
  }
}

int gyrotropic_susceptibility::num_cinternal_notowned_needed(component c,
                                                             void *P_internal_data) const {
  (void)c;
  (void)P_internal_data;
  return 0;
}

realnum *gyrotropic_susceptibility::cinternal_notowned_ptr(int inotowned, component c, int cmp,
                                                           int n, void *P_internal_data) const {
  gyrotropy_data *d = (gyrotropy_data *)P_internal_data;
  if (!d || !d->P[c][cmp][inotowned]) return NULL;
  return d->P[c][cmp][inotowned] + n;
}

void gyrotropic_susceptibility::dump_params(h5file *h5f, size_t *start) {
  size_t num_params = 9;
  size_t params_dims[1] = {num_params};
  realnum bias[] = {gyro_tensor[Y][Z], gyro_tensor[Z][X], gyro_tensor[X][Y]};
  realnum params_data[] = {8,     (realnum)get_id(), bias[X], bias[Y], bias[Z], omega_0, gamma,
                           alpha, (realnum)model};
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;
}

} // namespace meep
