# meep documentation map: Simulation Workflows

Use this map from `skills/meep-simulation-workflows/SKILL.md`.

## Workflow control docs
- `doc/docs/Python_User_Interface.md` - monitor APIs, run control, output routines.
- `doc/docs/Scheme_User_Interface.md` - Scheme run and monitor APIs.
- `doc/docs/The_Run_Function_Is_Not_A_Loop.md` - correct run-function usage.
- `doc/docs/Field_Functions.md` - step function and callback mechanics.
- `doc/docs/Synchronizing_the_Magnetic_and_Electric_Fields.md` - field synchronization details.
- `doc/docs/Parallel_Meep.md` - MPI execution and parallel behavior.
- `doc/docs/FAQ.md` - convergence checks and common diagnostics.

## Workflow tutorials
- `doc/docs/Python_Tutorials/Near_to_Far_Field_Spectra.md`
- `doc/docs/Scheme_Tutorials/Near_to_Far_Field_Spectra.md`
- `doc/docs/Python_Tutorials/Basics.md`

## Runnable workflow baselines
- `python/examples/bend-flux.py`
- `python/examples/cavity-farfield.py`
- `python/examples/antenna-radiation.py`
- `scheme/examples/bend-flux.ctl`
- `scheme/examples/cavity-farfield.ctl`
- `scheme/examples/antenna-radiation.ctl`

## Escalation
If workflow behavior remains unclear, use `skills/meep-simulation-workflows/references/source_map.md`.
