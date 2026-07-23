# gmat-mcp

⚠️ Note: This project is currently under development.

MCP (Model Context Protocol) server that lets Claude drive GMAT (General Mission
Analysis Tool, NASA) to build and analyze space missions.



## Requirements
- **GMAT** installed (tested with R2026a), with its shipped data: gravity
  coefficients (`data/gravity/earth/`), DE421 ephemerides
  (`data/planetary_ephem/`) and atmosphere models (`data/atmosphere/earth/`) —
  needed for spherical harmonics, third bodies and drag.
- **Python 3.10+** compatible with the GMAT binding (see `bin/gmatpy/` in your
  GMAT installation).
- **Python dependencies**: `pip install -r requirements.txt` in the venv
  (`mcp`, `matplotlib`, `plotly`, `numpy`, `sgp4`).
- **Local install files** (not versioned, machine-specific — see `.gitignore`):
  `load_gmat.py` and `api_startup_file.txt`, which point to the GMAT install
  path. Generate them locally during setup.

## Setup

Setup instructions coming soon.


## Structure
- `gmat_session.py` — wrapper around the GMAT API
- `server.py` — MCP server entry point
- `tools/` — MCP tool definitions

## Python dependencies
Installed in the `gmat-mcp-env` venv: `mcp`, `matplotlib`, `plotly`, `numpy`,
and `sgp4` (TLE import). Reinstall if needed:
`gmat-mcp-env\Scripts\python.exe -m pip install -r requirements.txt`.

## Available tools (33)
- **Session**: `set_mission`, `get_mission_info`, `reset_session`
- **Spacecraft**: `create_spacecraft` (Keplerian *or* Cartesian state),
  `get_spacecraft_state`, `list_spacecraft`,
  `create_spacecraft_from_tle`, `parse_tle`
- **Propagation**: `create_propagator` (point-mass by default; optional Jn
  harmonics, Moon/Sun third bodies, drag, SRP), `propagate`
- **Maneuvers**: `apply_impulsive_burn`, `create_thruster`,
  `apply_finite_burn`, `target_maneuver` (differential corrector)
- **Analysis / viz**: `plot_trajectory_3d`, `plot_orbital_elements`,
  `plot_ground_track`, `export_data`, `get_report`
- **Export**: `export_script` (GMAT `.script`), `export_python`
- **Analytic calculations**: `suggest_hohmann_transfer`,
  `suggest_bielliptic_transfer`, `combined_plane_change`, `solve_lambert`,
  `time_to_apsis`, `convert_time`, `keplerian_to_equinoctial`,
  `equinoctial_to_keplerian`, `tsiolkovsky`, `calculate_escape_velocity`,
  `validate_orbit`, `get_gravitational_parameters`

## Known limitations
The server drives GMAT through the direct API (`gmat.Construct`), without the
script/Sandbox engine. A few capabilities therefore rely on documented
approximations:

- **`apply_finite_burn`**: finite thrust is discretized into small impulses
  (symmetric splitting). Replaying the exported GMAT `.script` (a real
  `BeginFiniteBurn` sequence) may differ by ~0.1%.
- **`target_maneuver`**: home-grown differential corrector (shooting +
  bisection). It relies on a snapshot/restore of the propagator's Cartesian
  state; the converged burn is exported as an ordinary `ImpulsiveBurn`, but the
  `Target/Vary/Achieve` sequence is not regenerated in the `.script`. The
  restore does not roll back GMAT's internal epoch (negligible over a few
  iterations with near-stationary forces).
- **`create_spacecraft_from_tle`**: TEME→EarthMJ2000Eq conversion done by GMAT;
  the A1/UTC epoch is approximated (~37 s, negligible for the frame rotation).
  A TLE provides *mean* elements: the osculating SMA read back can differ by
  ~10 km, which is expected.
- **Drag / SRP**: constant solar and geomagnetic flux (F10.7 = 150), with no
  date dependence and no space-weather files.

## License
MIT — see [LICENSE](LICENSE).
