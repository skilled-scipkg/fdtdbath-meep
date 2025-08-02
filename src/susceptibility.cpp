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

#include <numeric>
#include <stdlib.h>
#include <string.h>
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
            pcg_extras::seed_seq_from<std::random_device> seed_source;
            pcg32_fast random(seed_source);
            cxx::ziggurat_normal_distribution<double> normal_distr(0, gaussian_random_amp);
            // I can also make this faster by generating num_bath gaussian random numbers at once
            /*
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
            */

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
