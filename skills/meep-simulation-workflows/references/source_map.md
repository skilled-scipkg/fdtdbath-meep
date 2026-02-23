# meep source map: Simulation Workflows

Use this map only after reading `skills/meep-simulation-workflows/references/doc_map.md`.

## Starter commands
- `python python/examples/bend-flux.py`
- `python python/examples/cavity-farfield.py`
- `meep scheme/examples/bend-flux.ctl`
- `mpirun -np 4 python -m mpi4py python/examples/bend-flux.py`

## Validation checkpoints
- Normalization and scattering runs use identical discretization/monitor placement.
- Decay-based termination (`stop_when_fields_decayed`) is satisfied before reading spectra.
- Flux/far-field outputs are finite and approximately conservation-consistent.
- MPI and serial runs agree within expected discretization/roundoff tolerances.

## Function-level source entry points
- `python/simulation.py` - `run`, `add_flux`, `add_near2far`, `get_dft_array`, `solve_cw`, `reset_meep`.
- `python/simulation.py` - step-function helpers `at_every`, `at_beginning`, `at_end`, `stop_when_fields_decayed`, `synchronized_magnetic`.
- `scheme/meep.scm.in` - `(run-until ...)`, `(run-sources ...)`, `(run-sources+ ...)`, `(at-every ...)`, `(synchronized-magnetic ...)`.
- `src/step.cpp` - `fields::step`, `fields::step_source`, `fields::step_boundaries`.
- `src/step_generic.cpp` - generic stepping path used in runtime loops.
- `src/step_db.cpp` - B/D update path (`fields::step_db`).
- `src/time.cpp` - wall-time and timing bookkeeping.
- `src/update_eh.cpp` - E/H update kernel dispatch.
- `src/dft.cpp` - flux/mode DFT accumulation and retrieval.
- `src/near2far.cpp` - near-to-far data accumulation and far-field evaluation.
- `src/h5fields.cpp` - field I/O behavior during workflow output.
- `src/loop_in_chunks.cpp` - chunked execution and parallel loop orchestration.
- `src/boundaries.cpp` - boundary conditions and Bloch handling.

## Fast source navigation
- `rg -n "def (run|add_flux|add_near2far|get_dft_array|solve_cw|reset_meep)|def (at_every|stop_when_fields_decayed|synchronized_magnetic)" python/simulation.py`
- `rg -n "define \(run-until|define \(run-sources|define \(run-sources\+|define \(at-every|define \(synchronized-magnetic" scheme/meep.scm.in`
- `rg -n "fields::step|step_db|update_eh|add_dft_flux|add_dft_near2far|loop_in_chunks|set_boundary" src`
