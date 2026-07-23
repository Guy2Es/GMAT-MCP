# gmat_session.py
"""
Wrapper around the GMAT API (load_gmat / gmatpy).
Centralizes the creation and tracking of GMAT objects (spacecraft,
propagators, etc.) created during the MCP server session.
"""
import load_gmat
from load_gmat import gmat
import math
import os

# Project base folder (where gmat_session.py itself lives), used to resolve
# relative export paths reliably, independently of the working directory
# (cwd) of the process that launches the MCP server (Claude Desktop does not
# guarantee a specific cwd).
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

EARTH_MU_KM3_S2 = 398600.4415  # Earth's standard gravitational parameter
EARTH_RADIUS_KM = 6378.137  # Earth's equatorial radius
EARTH_ROTATION_RAD_S = 7.2921158553e-5  # Earth's rotation rate (rad/s)

_MONTHS = {"Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
           "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12}


def _epoch_to_jd(epoch: str) -> float:
    """
    Converts a GMAT UTCGregorian epoch ('01 Jan 2026 00:00:00.000')
    to a Julian date. Manual parsing rather than strptime('%b') so as not
    to depend on the process locale (GMAT months always in English).
    """
    day_s, month_s, year_s, time_s = epoch.strip().split()
    day, month, year = int(day_s), _MONTHS[month_s], int(year_s)
    hh, mm, ss = time_s.split(":")
    day_frac = (int(hh) + int(mm) / 60 + float(ss) / 3600) / 24

    # Standard Gregorian -> JD algorithm (Meeus)
    if month <= 2:
        year -= 1
        month += 12
    a = year // 100
    b = 2 - a + a // 4
    jd0 = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.5
    return jd0 + day_frac


def _jd_to_gregorian(jd: float) -> str:
    """
    Converts a Julian date to a GMAT UTCGregorian string
    ('01 Jan 2026 00:00:00.000'). Inverse of _epoch_to_jd (Meeus algorithm).
    """
    jd += 0.5
    z = int(jd)
    f = jd - z
    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - alpha // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715

    day_frac = f * 24.0
    hh = int(day_frac)
    mm = int((day_frac - hh) * 60.0)
    ss = (((day_frac - hh) * 60.0) - mm) * 60.0
    month_name = [k for k, v in _MONTHS.items() if v == month][0]
    return f"{day:02d} {month_name} {year} {hh:02d}:{mm:02d}:{ss:06.3f}"


def _gmst_rad(jd_ut1: float) -> float:
    """
    Greenwich mean sidereal time (IAU 1982 formula, cf. Vallado).
    UTC is treated as UT1 (difference < 0.9 s, negligible here).
    """
    t = (jd_ut1 - 2451545.0) / 36525.0
    gmst_sec = (67310.54841
                + (876600.0 * 3600 + 8640184.812866) * t
                + 0.093104 * t**2
                - 6.2e-6 * t**3)
    return math.radians((gmst_sec % 86400.0) / 240.0)


def _cartesian_to_keplerian(state, mu=EARTH_MU_KM3_S2):
    """
    Converts a cartesian state [rx, ry, rz, vx, vy, vz] (km, km/s)
    to Keplerian orbital elements. Standard two-body formulas.
    """
    rx, ry, rz, vx, vy, vz = state
    r_vec = (rx, ry, rz)
    v_vec = (vx, vy, vz)

    r = math.sqrt(rx**2 + ry**2 + rz**2)
    v = math.sqrt(vx**2 + vy**2 + vz**2)

    # Angular momentum
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    h = math.sqrt(hx**2 + hy**2 + hz**2)

    # Line of nodes
    nx, ny, nz = -hy, hx, 0.0
    n = math.sqrt(nx**2 + ny**2)

    # Eccentricity vector
    rv_dot = rx * vx + ry * vy + rz * vz
    ex = (1 / mu) * ((v**2 - mu / r) * rx - rv_dot * vx)
    ey = (1 / mu) * ((v**2 - mu / r) * ry - rv_dot * vy)
    ez = (1 / mu) * ((v**2 - mu / r) * rz - rv_dot * vz)
    ecc = math.sqrt(ex**2 + ey**2 + ez**2)

    # Semi-major axis (via specific energy). For an escape orbit
    # (ecc >= 1, energy >= 0) the sma is negative (hyperbola) or infinite
    # (parabola): we leave it as-is rather than returning a misleading
    # positive value, and we classify the orbit type for the caller.
    energy = v**2 / 2 - mu / r
    if abs(energy) < 1e-12:
        sma = float("inf")  # parabola: semi-major axis undefined
        orbit_type = "parabolic"
    else:
        sma = -mu / (2 * energy)
        orbit_type = "hyperbolic" if ecc >= 1.0 else "elliptical"

    # Inclination
    inc = math.degrees(math.acos(max(-1.0, min(1.0, hz / h))))

    # RAAN
    if n > 1e-10:
        raan = math.degrees(math.acos(max(-1.0, min(1.0, nx / n))))
        if ny < 0:
            raan = 360 - raan
    else:
        raan = 0.0  # equatorial orbit, RAAN undefined

    # Argument of periapsis
    if n > 1e-10 and ecc > 1e-10:
        aop = math.degrees(math.acos(max(-1.0, min(1.0, (nx * ex + ny * ey) / (n * ecc)))))
        if ez < 0:
            aop = 360 - aop
    else:
        aop = 0.0

    # True anomaly
    if ecc > 1e-10:
        ta = math.degrees(math.acos(max(-1.0, min(1.0, (ex * rx + ey * ry + ez * rz) / (ecc * r)))))
        if rv_dot < 0:
            ta = 360 - ta
    else:
        ta = 0.0  # circular orbit, true anomaly classically undefined

    return {
        "sma_km": sma,
        "ecc": ecc,
        "inc_deg": inc,
        "raan_deg": raan,
        "aop_deg": aop,
        "ta_deg": ta,
        "orbit_type": orbit_type,
    }

class GmatSession:
    """
    Represents a single GMAT session, living for the entire lifetime
    of the server process.
    """

    def __init__(self):
        self.spacecraft = {}
        self.propagators = {}
        self.spacecraft_creation_state = {}
        self.thrusters = {}  # thrusters/tanks for FiniteBurn (Phase 2)
        self.mission_name = None
        self.output_dir = None  # None until set_mission() has been called
        # GMAT frames/converter built on demand for TLE import
        # (TEME -> EarthMJ2000Eq); reset to None by reset() (gmat.Clear).
        self._teme_cs = None
        self._eci_cs = None
        self._coord_conv = None

    def set_mission(self, name: str, base_path: str = None) -> dict:
        """
        Sets the current mission name and its output folder.
        All subsequent exports (export_script, export_python, export_data,
        plot_*) resolved via resolve_path() will go into this folder.
        """
        if not name or not name.strip():
            raise ValueError("The mission name cannot be empty.")

        # Sanitize: basename only, forbidden Windows characters removed --
        # the mission name feeds a folder path, we don't want it to be able
        # to escape it (path traversal) or break under Windows.
        sanitized = name.strip()
        for ch in '<>:"/\\|?*':
            sanitized = sanitized.replace(ch, "_")
        sanitized = os.path.basename(sanitized)
        if not sanitized or sanitized in (".", ".."):
            raise ValueError(f"Invalid mission name: '{name}'.")

        base = base_path or PROJECT_DIR
        output_dir = os.path.join(base, "output", sanitized)
        os.makedirs(output_dir, exist_ok=True)

        self.mission_name = sanitized
        self.output_dir = output_dir

        return {"mission_name": sanitized, "output_dir": output_dir}

    def resolve_path(self, filename: str) -> str:
        """
        Resolves an export file name to an absolute path.

        If a mission is active (set_mission() already called): only the file
        name (basename) is kept -- any path/folder provided by the caller is
        ignored -- and the file is placed in the current mission's folder,
        whatever was passed.

        Otherwise (legacy behavior, no mission defined): a relative path is
        resolved against PROJECT_DIR rather than the process cwd (not
        guaranteed under Claude Desktop); an absolute path is kept as-is.
        """
        if self.output_dir is not None:
            path = os.path.join(self.output_dir, os.path.basename(filename))
        elif os.path.isabs(filename):
            path = filename
        else:
            path = os.path.join(PROJECT_DIR, filename)

        os.makedirs(os.path.dirname(path), exist_ok=True)
        return path

    def reset(self) -> dict:
        """
        Resets the session: clears all GMAT objects (spacecraft,
        propagators, force models, burns...) and the Python tracking
        dictionaries, to start from a fresh session without having to
        restart the server. The mission name and output folder are KEPT
        (a mission reset is a distinct action, via set_mission()).

        Useful in particular to reuse an already-taken spacecraft/propagator
        name, or to chain several missions within the same process.
        """
        n_sats = len(self.spacecraft)
        n_props = len(self.propagators)

        # gmat.Clear() empties the engine-side GMAT object store; without it,
        # recreating an object of the same name would fail ("already exists").
        gmat.Clear()

        self.spacecraft = {}
        self.propagators = {}
        self.spacecraft_creation_state = {}
        self.thrusters = {}
        self._teme_cs = None
        self._eci_cs = None
        self._coord_conv = None

        return {
            "reset": True,
            "spacecraft_cleared": n_sats,
            "propagators_cleared": n_props,
            "mission_name": self.mission_name,
            "output_dir": self.output_dir,
        }

    def create_spacecraft(self, name: str, sma_km: float = None, ecc: float = None,
                           inc_deg: float = None, raan_deg: float = 0.0,
                           aop_deg: float = 0.0, ta_deg: float = 0.0,
                           epoch: str = "01 Jan 2026 00:00:00.000",
                           cartesian_state: list = None) -> dict:
        """
        Creates a spacecraft, either from Keplerian elements (sma_km/ecc/inc_deg...)
        or from a cartesian state [x, y, z, vx, vy, vz] (km, km/s) in the
        inertial EarthMJ2000Eq frame. The two modes are mutually exclusive.
        """
        if name in self.spacecraft:
            raise ValueError(
                f"A spacecraft named '{name}' already exists in this session."
            )

        use_cartesian = cartesian_state is not None
        if use_cartesian:
            if sma_km is not None or ecc is not None or inc_deg is not None:
                raise ValueError(
                    "cartesian_state and the Keplerian elements (sma_km/ecc/inc_deg) "
                    "are mutually exclusive."
                )
            if len(cartesian_state) != 6:
                raise ValueError(
                    "cartesian_state must contain 6 values [x, y, z, vx, vy, vz]."
                )
        elif sma_km is None or ecc is None or inc_deg is None:
            raise ValueError(
                "Provide either the Keplerian elements (sma_km, ecc, inc_deg), "
                "or a cartesian_state."
            )

        sat = gmat.Construct("Spacecraft", name)
        sat.SetField("DateFormat", "UTCGregorian")
        sat.SetField("Epoch", epoch)
        sat.SetField("CoordinateSystem", "EarthMJ2000Eq")

        if use_cartesian:
            x, y, z, vx, vy, vz = cartesian_state
            sat.SetField("DisplayStateType", "Cartesian")
            sat.SetField("X", x)
            sat.SetField("Y", y)
            sat.SetField("Z", z)
            sat.SetField("VX", vx)
            sat.SetField("VY", vy)
            sat.SetField("VZ", vz)
            # Creation state stored in Keplerian form (the .script/.py export
            # declares spacecraft in Keplerian, cf. export_script).
            kep = _cartesian_to_keplerian(cartesian_state)
            creation = {
                "sma_km": kep["sma_km"],
                "ecc": kep["ecc"],
                "inc_deg": kep["inc_deg"],
                "raan_deg": kep["raan_deg"],
                "aop_deg": kep["aop_deg"],
                "ta_deg": kep["ta_deg"],
                "epoch": epoch,
            }
        else:
            sat.SetField("DisplayStateType", "Keplerian")
            sat.SetField("SMA", sma_km)
            sat.SetField("ECC", ecc)
            sat.SetField("INC", inc_deg)
            sat.SetField("RAAN", raan_deg)
            sat.SetField("AOP", aop_deg)
            sat.SetField("TA", ta_deg)
            creation = {
                "sma_km": sma_km,
                "ecc": ecc,
                "inc_deg": inc_deg,
                "raan_deg": raan_deg,
                "aop_deg": aop_deg,
                "ta_deg": ta_deg,
                "epoch": epoch,
            }

        self.spacecraft[name] = sat

        # Stored separately from last_keplerian_state (Phase 2/3): export_script()
        # must declare the spacecraft in its CREATION state, not its final state
        # -- the Propagate/Maneuver timeline (Phase 5) is what is meant to
        # evolve the state from that starting point. Reusing last_keplerian_state
        # for the initial declaration would replay the entire mission on top of
        # an already-arrived spacecraft, doubling the effect of maneuvers/propagations.
        self.spacecraft_creation_state[name] = creation

        # Necessary for GMAT to convert the input state into an internal state
        # -- without this, a re-read (GetField) returns 0.
        gmat.Initialize()

        return {"name": name, **creation}

    @staticmethod
    def parse_tle(line1: str, line2: str) -> dict:
        """
        Reads a TLE (2 lines) and returns its mean elements (SGP4), without
        creating anything in GMAT. The semi-major axis is derived from the
        mean motion.
        """
        from sgp4.api import Satrec
        sat = Satrec.twoline2rv(line1.strip(), line2.strip())

        n_rad_min = sat.no_kozai      # mean motion (rad/min)
        n_rad_s = n_rad_min / 60.0
        sma_km = (EARTH_MU_KM3_S2 / n_rad_s**2) ** (1.0 / 3.0)
        period_sec = 2 * math.pi / n_rad_s

        jd_epoch = sat.jdsatepoch + sat.jdsatepochF
        return {
            "norad_id": sat.satnum,
            "epoch_utc": _jd_to_gregorian(jd_epoch),
            "epoch_jd": jd_epoch,
            "inc_deg": math.degrees(sat.inclo),
            "raan_deg": math.degrees(sat.nodeo),
            "ecc": sat.ecco,
            "aop_deg": math.degrees(sat.argpo),
            "mean_anomaly_deg": math.degrees(sat.mo),
            "mean_motion_rev_day": n_rad_min * 1440.0 / (2 * math.pi),
            "sma_km": sma_km,
            "period_sec": period_sec,
            "bstar": sat.bstar,
        }

    def _ensure_teme_frames(self):
        """Builds (once) the TEME/MJ2000Eq frames and the converter."""
        if self._coord_conv is not None:
            return
        ss = gmat.GetSolarSystem()
        self._teme_cs = gmat.Construct("CoordinateSystem", "TLE_TEME", "Earth", "TEME")
        self._eci_cs = gmat.Construct("CoordinateSystem", "TLE_MJ2000", "Earth", "MJ2000Eq")
        self._teme_cs.SetSolarSystem(ss)
        self._eci_cs.SetSolarSystem(ss)
        gmat.Initialize()
        self._coord_conv = gmat.CoordinateConverter()

    def create_spacecraft_from_tle(self, name: str, line1: str, line2: str) -> dict:
        """
        Creates a spacecraft from a TLE: propagates SGP4 to the TLE epoch to
        obtain the state (position/velocity) in the TEME frame, converts it to
        EarthMJ2000Eq via GMAT, and creates the spacecraft at that cartesian state.

        Call BEFORE create_propagator (like create_spacecraft, the method
        calls gmat.Initialize()). The osculating/mean difference means the
        re-read SMA may differ by ~10 km from the TLE's mean elements — this is normal.
        """
        from sgp4.api import Satrec
        sat = Satrec.twoline2rv(line1.strip(), line2.strip())

        jd, fr = sat.jdsatepoch, sat.jdsatepochF
        err, r_teme, v_teme = sat.sgp4(jd, fr)
        if err != 0:
            raise ValueError(f"SGP4 propagation failed (code {err}) at the TLE epoch.")

        self._ensure_teme_frames()

        # A1ModJulian epoch ≈ JD_utc - 2430000.0. The A1/UTC difference (~37 s)
        # is negligible for the TEME->J2000 rotation matrix (precession/
        # nutation vary on the order of arcsec/day).
        epoch_a1mjd = (jd + fr) - 2430000.0
        instate = gmat.Rvector6(r_teme[0], r_teme[1], r_teme[2],
                                v_teme[0], v_teme[1], v_teme[2])
        outstate = gmat.Rvector6()
        self._coord_conv.Convert(epoch_a1mjd, instate, self._teme_cs,
                                 outstate, self._eci_cs)
        j2000_state = [outstate[i] for i in range(6)]

        epoch_str = _jd_to_gregorian(jd + fr)
        result = self.create_spacecraft(name, epoch=epoch_str,
                                        cartesian_state=j2000_state)
        result["norad_id"] = sat.satnum
        result["source"] = "TLE/SGP4"
        return result

    def get_spacecraft_state(self, name: str) -> dict:
        if name not in self.spacecraft:
            raise ValueError(f"No spacecraft named '{name}' in this session.")

        sat = self.spacecraft[name]
        return {
            "name": name,
            "sma_km": float(sat.GetField("SMA")),
            "ecc": float(sat.GetField("ECC")),
            "inc_deg": float(sat.GetField("INC")),
            "raan_deg": float(sat.GetField("RAAN")),
            "aop_deg": float(sat.GetField("AOP")),
            "ta_deg": float(sat.GetField("TA")),
        }

    def list_spacecraft(self) -> list:
        return list(self.spacecraft.keys())

    # Gravitational coefficient files available in the GMAT installation
    # (C:/GMAT/data/gravity/earth/) and their maximum degree/order.
    _GRAVITY_FILES = {
        "JGM2": ("JGM2.cof", 70),
        "JGM3": ("JGM3.cof", 70),
        "EGM96": ("EGM96.cof", 360),
    }
    # Bodies accepted as third-bodies (tolerant of the usual LLM aliases).
    _THIRD_BODY_ALIASES = {
        "MOON": "Luna", "LUNA": "Luna",
        "SUN": "Sun", "SOL": "Sun",
        "JUPITER": "Jupiter", "MARS": "Mars", "VENUS": "Venus",
    }

    def create_propagator(self, spacecraft_name: str,
                           propagator_name: str = None,
                           step_size_sec: float = 60.0,
                           gravity_degree: int = 0,
                           gravity_order: int = None,
                           gravity_file: str = "JGM3",
                           third_bodies: list = None,
                           enable_drag: bool = False,
                           atmosphere_model: str = "MSISE90",
                           enable_srp: bool = False,
                           dry_mass_kg: float = 850.0,
                           cd: float = 2.2,
                           drag_area_m2: float = 15.0,
                           cr: float = 1.8,
                           srp_area_m2: float = 1.0) -> dict:
        """
        Creates a ForceModel + Propagator associated with an existing spacecraft.

        By default (gravity_degree=0), the force model is a simple Earth
        point-mass — identical to the historical behavior. Pass parameters to
        enrich the model (harmonics, third-bodies, drag, SRP).
        """
        if spacecraft_name not in self.spacecraft:
            raise ValueError(f"No spacecraft named '{spacecraft_name}' in this session.")

        prop_name = propagator_name or f"Prop_{spacecraft_name}"
        if prop_name in self.propagators:
            raise ValueError(
                f"A propagator named '{prop_name}' already exists in this session."
            )

        sat = self.spacecraft[spacecraft_name]

        gravity_order = gravity_degree if gravity_order is None else gravity_order

        # --- Force-parameter validation (before any construction) ---
        third_body_names = []
        if third_bodies:
            for body in third_bodies:
                key = str(body).strip().upper()
                if key not in self._THIRD_BODY_ALIASES:
                    raise ValueError(
                        f"Third-body '{body}' unknown — expected one of "
                        f"{sorted(set(self._THIRD_BODY_ALIASES.values()))}."
                    )
                third_body_names.append(self._THIRD_BODY_ALIASES[key])

        if gravity_degree > 0:
            if gravity_file not in self._GRAVITY_FILES:
                raise ValueError(
                    f"Gravity file '{gravity_file}' unknown — expected one of "
                    f"{sorted(self._GRAVITY_FILES)}."
                )
            _, max_degree = self._GRAVITY_FILES[gravity_file]
            if gravity_degree > max_degree:
                raise ValueError(
                    f"gravity_degree={gravity_degree} exceeds the maximum "
                    f"({max_degree}) of the {gravity_file} model."
                )
            if gravity_order > gravity_degree:
                raise ValueError(
                    f"gravity_order={gravity_order} cannot exceed "
                    f"gravity_degree={gravity_degree}."
                )

        # --- Spacecraft physical properties (only if drag/SRP) ---
        # Must be set BEFORE gmat.Initialize() / PrepareInternals().
        if enable_drag:
            sat.SetField("DryMass", dry_mass_kg)
            sat.SetField("Cd", cd)
            sat.SetField("DragArea", drag_area_m2)
        if enable_srp:
            sat.SetField("DryMass", dry_mass_kg)
            sat.SetField("Cr", cr)
            sat.SetField("SRPArea", srp_area_m2)

        # --- Force model ---
        fm = gmat.Construct("ForceModel", f"FM_{spacecraft_name}")
        fm.SetField("CentralBody", "Earth")

        # Earth gravity: spherical harmonics if a degree is requested,
        # otherwise point-mass (historical behavior).
        if gravity_degree > 0:
            cof_file, _ = self._GRAVITY_FILES[gravity_file]
            earthgrav = gmat.Construct("GravityField")
            earthgrav.SetField("BodyName", "Earth")
            earthgrav.SetField("Degree", gravity_degree)
            earthgrav.SetField("Order", gravity_order)
            earthgrav.SetField("PotentialFile", cof_file)
            fm.AddForce(earthgrav)
        else:
            pointmass = gmat.Construct("PointMassForce")
            pointmass.SetField("BodyName", "Earth")
            fm.AddForce(pointmass)

        # Third-bodies: one PointMassForce per body (DE421 ephemerides).
        for body in third_body_names:
            tb = gmat.Construct("PointMassForce")
            tb.SetField("BodyName", body)
            fm.AddForce(tb)

        # Atmospheric drag.
        if enable_drag:
            drag = gmat.Construct("DragForce")
            drag.SetField("AtmosphereModel", atmosphere_model)
            atmos = gmat.Construct(atmosphere_model)
            drag.SetReference(atmos)
            # Constant solar/geomagnetic flux (v1): avoids any dependence on
            # the date and on space-weather files.
            drag.SetField("HistoricWeatherSource", "ConstantFluxAndGeoMag")
            drag.SetField("PredictedWeatherSource", "ConstantFluxAndGeoMag")
            drag.SetField("F107", 150.0)
            drag.SetField("F107A", 150.0)
            drag.SetField("MagneticIndex", 3.0)
            fm.AddForce(drag)

        # Solar radiation pressure.
        if enable_srp:
            srp = gmat.Construct("SolarRadiationPressure")
            srp.SetField("SRPModel", "Spherical")
            fm.AddForce(srp)

        # --- Propagator + integrator ---
        pdprop = gmat.Construct("Propagator", prop_name)
        integrator = gmat.Construct("PrinceDormand78", f"Integ_{spacecraft_name}")
        pdprop.SetReference(integrator)
        pdprop.SetReference(fm)
        pdprop.SetField("InitialStepSize", step_size_sec)
        pdprop.SetField("Accuracy", 1.0e-12)
        pdprop.SetField("MinStep", 0.0)

        gmat.Initialize()

        pdprop.AddPropObject(sat)
        pdprop.PrepareInternals()

        # Important: the real object to "step" is not pdprop itself, but this
        # reference retrieved after PrepareInternals().
        gator = pdprop.GetPropagator()

        # Configuration stored for traceability and the .script/.py export.
        force_config = {
            "gravity_degree": gravity_degree,
            "gravity_order": gravity_order,
            "gravity_file": gravity_file,
            "third_bodies": third_body_names,
            "enable_drag": enable_drag,
            "atmosphere_model": atmosphere_model,
            "enable_srp": enable_srp,
            "dry_mass_kg": dry_mass_kg,
            "cd": cd,
            "drag_area_m2": drag_area_m2,
            "cr": cr,
            "srp_area_m2": srp_area_m2,
        }

        self.propagators[prop_name] = {
            "pdprop": pdprop,
            "gator": gator,
            "step_size_sec": step_size_sec,
            "spacecraft_name": spacecraft_name,
            "force_config": force_config,
        }

        return {
            "propagator_name": prop_name,
            "spacecraft_name": spacecraft_name,
            "step_size_sec": step_size_sec,
            "force_model": force_config,
        }

    def propagate(self, propagator_name: str, duration_sec: float) -> dict:
        """
        Propagates the orbit associated with this propagator for duration_sec
        seconds, and returns the final orbital state (converted from the raw
        cartesian vector, see note on resynchronization).
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")

        entry = self.propagators[propagator_name]
        gator = entry["gator"]
        step = entry["step_size_sec"]

        # History sampled at each Step(), on a global time axis (accumulated
        # over all propagate() calls of this propagator) — it is the raw
        # material of the visualization tools (Phase 6).
        base_elapsed = entry.get("total_propagated_sec", 0.0)
        history = entry.setdefault("history", [])
        if not history:
            history.append({
                "elapsed_sec": base_elapsed,
                "cartesian": list(gator.GetState()),
            })

        elapsed = 0.0
        while elapsed < duration_sec:
            this_step = min(step, duration_sec - elapsed)
            gator.Step(this_step)
            elapsed += this_step
            history.append({
                "elapsed_sec": base_elapsed + elapsed,
                "cartesian": list(gator.GetState()),
            })

        final_cartesian_state = gator.GetState()
        keplerian = _cartesian_to_keplerian(final_cartesian_state)

        # NEW: we store the result so that get_report()/export_script()
        # reflect the actually-propagated state, not just the creation state.
        entry["total_propagated_sec"] = entry.get("total_propagated_sec", 0.0) + duration_sec
        entry["last_keplerian_state"] = keplerian
        entry["last_cartesian_state"] = list(final_cartesian_state)
        entry["last_state_source"] = "propagated"
        entry.setdefault("timeline", []).append(
            {"type": "propagate", "duration_sec": duration_sec}
        )

        return {
            "name": entry["spacecraft_name"],
            **keplerian,
        }

    def apply_impulsive_burn(self, propagator_name: str, delta_v1_km_s: float,
                              delta_v2_km_s: float = 0.0, delta_v3_km_s: float = 0.0,
                              frame: str = "VNB") -> dict:
        """
        Applies an impulsive maneuver (instantaneous delta-v) to the
        spacecraft associated with this propagator, in the requested frame,
        and returns the resulting orbital state.

        Args:
            frame: "VNB" (prograde/normal/binormal, default), "LVLH",
                "MJ2000Eq" (these three in an Earth-centered local frame), or
                "Inertial" (delta-v expressed directly in EarthMJ2000Eq).
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")

        local_axes = {"VNB", "LVLH", "MJ2000Eq"}
        if frame not in local_axes and frame != "Inertial":
            raise ValueError(f"Frame '{frame}' unknown — expected one of {sorted(local_axes)} or 'Inertial'.")

        entry = self.propagators[propagator_name]
        sat_name = entry["spacecraft_name"]
        sat = self.spacecraft[sat_name]
        gator = entry["gator"]

        # PITFALL 1 (found in practice): sat has a cartesian state distinct
        # from gator's (same desynchronization as after Step(), cf. Phase 2)
        # -- without this sync, Fire() applies the delta-v to the CREATION
        # state, silently, even if the spacecraft has already been propagated.
        gator.UpdateSpaceObject()

        burn_count = sum(1 for t in entry.get("timeline", []) if t["type"] == "maneuver")
        burn_name = f"Burn_{propagator_name}_{burn_count + 1}"
        burn = gmat.Construct("ImpulsiveBurn", burn_name)

        if frame == "Inertial":
            burn.SetField("CoordinateSystem", "EarthMJ2000Eq")
        else:
            burn.SetField("CoordinateSystem", "Local")
            burn.SetField("Origin", "Earth")
            burn.SetField("Axes", frame)

        burn.SetField("Element1", delta_v1_km_s)
        burn.SetField("Element2", delta_v2_km_s)
        burn.SetField("Element3", delta_v3_km_s)

        # PITFALL 2 (found in practice): without SetSolarSystem, Fire() raises
        # "SolarSystem is NULL". Do NOT fix this with a second
        # gmat.Initialize() -- it invalidates the `gator` handle already
        # obtained via GetPropagator() (confirmed: GetState() then raises a
        # RuntimeError). SetSolarSystem is the targeted fix, with no side
        # effect on the rest of the session.
        burn.SetSolarSystem(gmat.GetSolarSystem())
        burn.SetSpacecraftToManeuver(sat)
        burn.Fire()

        # Push sat's post-burn state back to gator, so that the next Step()
        # calls start from the correct state.
        gator.UpdateFromSpaceObject()

        # Post-burn sample in the history (same elapsed_sec: the burn is
        # instantaneous) — the velocity discontinuity thus appears in the
        # plots at the right place.
        at_elapsed_sec = entry.get("total_propagated_sec", 0.0)
        entry.setdefault("history", []).append({
            "elapsed_sec": at_elapsed_sec,
            "cartesian": list(gator.GetState()),
        })

        keplerian = {
            "sma_km": float(sat.GetField("SMA")),
            "ecc": float(sat.GetField("ECC")),
            "inc_deg": float(sat.GetField("INC")),
            "raan_deg": float(sat.GetField("RAAN")),
            "aop_deg": float(sat.GetField("AOP")),
            "ta_deg": float(sat.GetField("TA")),
        }

        entry["last_keplerian_state"] = keplerian
        entry["last_state_source"] = "maneuvered"
        entry.setdefault("timeline", []).append({
            "type": "maneuver",
            "burn_name": burn_name,
            "frame": frame,
            "delta_v_km_s": [delta_v1_km_s, delta_v2_km_s, delta_v3_km_s],
            "at_elapsed_sec": at_elapsed_sec,
        })

        return {
            "name": sat_name,
            "burn_name": burn_name,
            "frame": frame,
            "delta_v_km_s": [delta_v1_km_s, delta_v2_km_s, delta_v3_km_s],
            **keplerian,
        }

    def create_thruster(self, spacecraft_name: str, thruster_name: str = None,
                        thrust_n: float = 440.0, isp_s: float = 320.0,
                        fuel_mass_kg: float = 200.0,
                        tank_name: str = None) -> dict:
        """
        Declares a chemical thruster (constant thrust) + its tank for a
        spacecraft, in view of a finite-burn maneuver (apply_finite_burn).

        The configuration is stored on the Python side (thrust, Isp,
        propellant): the burn dynamics are integrated by us segment by
        segment. The GMAT ChemicalTank/ChemicalThruster objects are NOT
        constructed directly here, to avoid risking invalidating an
        already-prepared propagator (any later gmat.Initialize() breaks the
        internal handle); they are, however, faithfully generated in export_script().
        """
        if spacecraft_name not in self.spacecraft:
            raise ValueError(f"No spacecraft named '{spacecraft_name}' in this session.")

        thr_name = thruster_name or f"Thruster_{spacecraft_name}"
        if thr_name in self.thrusters:
            raise ValueError(f"A thruster named '{thr_name}' already exists in this session.")

        tnk_name = tank_name or f"Tank_{spacecraft_name}"

        self.thrusters[thr_name] = {
            "spacecraft_name": spacecraft_name,
            "tank_name": tnk_name,
            "thrust_n": thrust_n,
            "isp_s": isp_s,
            "fuel_mass_kg": fuel_mass_kg,
            "fuel_remaining_kg": fuel_mass_kg,
        }

        return {
            "thruster_name": thr_name,
            "tank_name": tnk_name,
            "spacecraft_name": spacecraft_name,
            "thrust_n": thrust_n,
            "isp_s": isp_s,
            "fuel_mass_kg": fuel_mass_kg,
        }

    # Thrust direction expressed in the local VNB frame
    # (Velocity/Normal/Binormal) -> (Element1, Element2, Element3).
    _FINITE_BURN_DIRECTIONS = {
        "prograde": (1.0, 0.0, 0.0),
        "retrograde": (-1.0, 0.0, 0.0),
        "normal": (0.0, 1.0, 0.0),
        "antinormal": (0.0, -1.0, 0.0),
    }

    def apply_finite_burn(self, propagator_name: str, duration_sec: float,
                          direction: str = "prograde",
                          thruster_name: str = None,
                          segment_sec: float = None) -> dict:
        """
        Applies a finite-burn maneuver (continuous thrust) for duration_sec,
        discretizing the thrust into small impulses.

        Approach: the duration is split into segments; on each segment we
        propagate half the step, apply the equivalent impulse
        (dv = thrust/mass · dt) via the proven ImpulsiveBurn mechanism, then
        propagate the other half (symmetric splitting, order 2). Mass and
        propellant consumption are tracked on the Python side (mdot = F/(Isp·g0)).

        Args:
            direction: "prograde" (default), "retrograde", "normal" or
                "antinormal" (local VNB frame).
            thruster_name: thruster to use (default: the spacecraft's only
                thruster, created via create_thruster).
            segment_sec: duration of a segment (default: ~ duration/20, capped
                by the propagator step).
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")
        if duration_sec <= 0:
            raise ValueError("duration_sec must be strictly positive.")
        if direction not in self._FINITE_BURN_DIRECTIONS:
            raise ValueError(
                f"Direction '{direction}' unknown — expected one of "
                f"{sorted(self._FINITE_BURN_DIRECTIONS)}."
            )

        entry = self.propagators[propagator_name]
        sat_name = entry["spacecraft_name"]
        sat = self.spacecraft[sat_name]
        gator = entry["gator"]

        # Thruster selection.
        if thruster_name is None:
            candidates = [n for n, t in self.thrusters.items()
                          if t["spacecraft_name"] == sat_name]
            if not candidates:
                raise ValueError(
                    f"No thruster for '{sat_name}' — call create_thruster() first."
                )
            if len(candidates) > 1:
                raise ValueError(
                    f"Several thrusters for '{sat_name}' — specify thruster_name "
                    f"among {candidates}."
                )
            thruster_name = candidates[0]
        elif thruster_name not in self.thrusters:
            raise ValueError(f"No thruster named '{thruster_name}' in this session.")

        thr = self.thrusters[thruster_name]
        if thr["spacecraft_name"] != sat_name:
            raise ValueError(
                f"The thruster '{thruster_name}' belongs to '{thr['spacecraft_name']}', "
                f"not to '{sat_name}'."
            )

        F = thr["thrust_n"]
        isp = thr["isp_s"]
        g0 = 9.80665
        mdot = F / (isp * g0)

        try:
            dry_mass = float(sat.GetField("DryMass"))
        except Exception:
            dry_mass = 850.0
        fuel = thr["fuel_remaining_kg"]
        if fuel <= 0:
            raise ValueError(f"The tank of thruster '{thruster_name}' is empty.")

        # Number of segments.
        default_seg = min(entry["step_size_sec"], duration_sec / 20.0)
        seg = segment_sec or default_seg
        seg = max(seg, 1e-6)
        n_seg = max(1, int(math.ceil(duration_sec / seg)))
        seg = duration_sec / n_seg

        ux, uy, uz = self._FINITE_BURN_DIRECTIONS[direction]

        # A single ImpulsiveBurn reused for all segments.
        fb_count = sum(1 for t in entry.get("timeline", []) if t["type"] == "finite_burn")
        burn_name = f"FiniteBurn_{propagator_name}_{fb_count + 1}"
        burn = gmat.Construct("ImpulsiveBurn", f"_seg_{burn_name}")
        burn.SetField("CoordinateSystem", "Local")
        burn.SetField("Origin", "Earth")
        burn.SetField("Axes", "VNB")
        burn.SetSolarSystem(gmat.GetSolarSystem())
        burn.SetSpacecraftToManeuver(sat)

        base_elapsed = entry.get("total_propagated_sec", 0.0)
        history = entry.setdefault("history", [])
        if not history:
            history.append({
                "elapsed_sec": base_elapsed,
                "cartesian": list(gator.GetState()),
            })

        dv_total = 0.0
        fuel_used = 0.0
        elapsed_local = 0.0
        fuel_exhausted = False

        for _ in range(n_seg):
            # Symmetric splitting: half-step / impulse / half-step.
            gator.Step(seg / 2.0)

            dm = mdot * seg
            if dm >= fuel:
                dm = fuel
                fuel_exhausted = True
            m_mid = dry_mass + fuel - dm / 2.0
            dv_seg = (F / m_mid) * seg / 1000.0  # km/s

            gator.UpdateSpaceObject()
            burn.SetField("Element1", ux * dv_seg)
            burn.SetField("Element2", uy * dv_seg)
            burn.SetField("Element3", uz * dv_seg)
            burn.Fire()
            gator.UpdateFromSpaceObject()

            gator.Step(seg / 2.0)

            fuel -= dm
            fuel_used += dm
            dv_total += dv_seg
            elapsed_local += seg
            history.append({
                "elapsed_sec": base_elapsed + elapsed_local,
                "cartesian": list(gator.GetState()),
            })
            if fuel_exhausted:
                break

        thr["fuel_remaining_kg"] = fuel

        keplerian = _cartesian_to_keplerian(gator.GetState())
        entry["total_propagated_sec"] = base_elapsed + elapsed_local
        entry["last_keplerian_state"] = keplerian
        entry["last_cartesian_state"] = list(gator.GetState())
        entry["last_state_source"] = "finite_burn"
        entry.setdefault("timeline", []).append({
            "type": "finite_burn",
            "burn_name": burn_name,
            "thruster_name": thruster_name,
            "tank_name": thr["tank_name"],
            "direction": direction,
            "direction_vec": [ux, uy, uz],
            "duration_sec": elapsed_local,
            "thrust_n": F,
            "isp_s": isp,
            "at_elapsed_sec": base_elapsed,
        })

        result = {
            "name": sat_name,
            "burn_name": burn_name,
            "thruster_name": thruster_name,
            "direction": direction,
            "duration_sec": elapsed_local,
            "delta_v_total_km_s": dv_total,
            "fuel_used_kg": fuel_used,
            "fuel_remaining_kg": fuel,
            "final_mass_kg": dry_mass + fuel,
            "n_segments": n_seg,
            **keplerian,
        }
        if fuel_exhausted:
            result["warning"] = (
                "Tank exhausted before the end of the requested duration — "
                "burn interrupted."
            )
        return result

    def _snapshot_state(self, entry: dict) -> dict:
        """
        Captures the current state of a propagator (gator's cartesian vector,
        elapsed time) so it can be returned to after a targeting shot.
        """
        gator = entry["gator"]
        return {
            "cartesian": list(gator.GetState()),
            "total_propagated_sec": entry.get("total_propagated_sec", 0.0),
        }

    def _restore_state(self, entry: dict, snap: dict) -> None:
        """
        Restores a state captured by _snapshot_state: rewrites the cartesian
        vector into the spacecraft then pushes it back to the propagator.
        Uses the same mechanism as the post-burn (UpdateFromSpaceObject) —
        crucially NOT gmat.Initialize(), which would invalidate the gator handle.

        Assumed limitation: the internal GMAT epoch is not rolled back; with no
        notable effect over a few targeting iterations for quasi-stationary
        forces.
        """
        sat = self.spacecraft[entry["spacecraft_name"]]
        x, y, z, vx, vy, vz = snap["cartesian"]
        sat.SetField("DisplayStateType", "Cartesian")
        sat.SetField("X", x)
        sat.SetField("Y", y)
        sat.SetField("Z", z)
        sat.SetField("VX", vx)
        sat.SetField("VY", vy)
        sat.SetField("VZ", vz)
        entry["gator"].UpdateFromSpaceObject()
        entry["total_propagated_sec"] = snap["total_propagated_sec"]

    _TARGET_GOALS = {
        "apoapsis_radius_km", "periapsis_radius_km",
        "apoapsis_altitude_km", "periapsis_altitude_km",
        "sma_km", "ecc", "period_sec",
    }
    _TARGET_CONTROL_AXIS = {"prograde": 1, "normal": 2, "binormal": 3}

    def target_maneuver(self, propagator_name: str, goal_type: str,
                        goal_value: float, control: str = "prograde",
                        initial_guess_km_s: float = 0.1,
                        tolerance: float = 1e-3,
                        max_iterations: int = 60) -> dict:
        """
        In-house differential corrector (shooting + secant method): searches
        for the delta-v that drives a final orbital quantity to a target value,
        then applies the converged maneuver.

        The targeting iterates "dry" (snapshot/restore of the state, no effect
        on the history or the timeline); only the converged burn is actually
        applied at the end, like an ordinary ImpulsiveBurn.

        Args:
            goal_type: target quantity — "apoapsis_radius_km",
                "periapsis_radius_km", "apoapsis_altitude_km",
                "periapsis_altitude_km", "sma_km", "ecc" or "period_sec".
            goal_value: target value (km, dimensionless or s depending on goal_type).
            control: delta-v control axis — "prograde" (default),
                "normal" or "binormal" (VNB frame).
            initial_guess_km_s: first delta-v trial (km/s).
            tolerance: tolerance on the residual (same unit as goal_value).
            max_iterations: maximum number of iterations.
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")
        if goal_type not in self._TARGET_GOALS:
            raise ValueError(
                f"goal_type '{goal_type}' unknown — expected one of {sorted(self._TARGET_GOALS)}."
            )
        if control not in self._TARGET_CONTROL_AXIS:
            raise ValueError(
                f"control '{control}' unknown — expected one of "
                f"{sorted(self._TARGET_CONTROL_AXIS)}."
            )

        entry = self.propagators[propagator_name]
        sat_name = entry["spacecraft_name"]
        sat = self.spacecraft[sat_name]
        gator = entry["gator"]
        mu = EARTH_MU_KM3_S2
        axis = self._TARGET_CONTROL_AXIS[control]

        # Finite sentinel (rather than an inf) for quantities undefined on an
        # escape orbit: apoapsis and period "go to infinity" as the escape
        # velocity is approached. Returning a large finite value cleanly backs
        # off the root finder instead of breaking its arithmetic.
        BIG = 1e15

        def measure(kep):
            sma = kep["sma_km"]
            ecc = kep["ecc"]
            if goal_type == "sma_km":
                return sma
            if goal_type == "ecc":
                return ecc
            if goal_type == "period_sec":
                if sma <= 0:
                    return BIG  # escape orbit: period undefined
                return 2 * math.pi * math.sqrt(sma**3 / mu)
            ra = sma * (1 + ecc)
            rp = sma * (1 - ecc)
            if goal_type == "apoapsis_radius_km":
                return BIG if sma <= 0 else ra
            if goal_type == "periapsis_radius_km":
                return rp
            if goal_type == "apoapsis_altitude_km":
                return BIG if sma <= 0 else ra - EARTH_RADIUS_KM
            return rp - EARTH_RADIUS_KM  # periapsis_altitude_km

        # Snapshot of the starting state (the gator must be synchronized).
        gator.UpdateSpaceObject()
        snap = self._snapshot_state(entry)

        burn = gmat.Construct("ImpulsiveBurn", f"_target_{propagator_name}")
        burn.SetField("CoordinateSystem", "Local")
        burn.SetField("Origin", "Earth")
        burn.SetField("Axes", "VNB")
        burn.SetSolarSystem(gmat.GetSolarSystem())
        burn.SetSpacecraftToManeuver(sat)

        def trial(dv):
            self._restore_state(entry, snap)
            gator.UpdateSpaceObject()
            burn.SetField("Element1", dv if axis == 1 else 0.0)
            burn.SetField("Element2", dv if axis == 2 else 0.0)
            burn.SetField("Element3", dv if axis == 3 else 0.0)
            burn.Fire()
            gator.UpdateFromSpaceObject()
            kep = _cartesian_to_keplerian(gator.GetState())
            return measure(kep) - goal_value

        # --- Bracketing a root, then false position (Illinois) ---
        # The targeted quantity is generally monotonic in the delta-v over the
        # bound domain; we look for an interval [a, b] where the residual
        # changes sign, starting from dv=0 and expanding from initial_guess.
        a = 0.0
        fa = trial(a)
        best_res = abs(fa)
        best_dv = a
        converged = False
        iterations = 0

        b = initial_guess_km_s if abs(initial_guess_km_s) > 1e-12 else 0.1
        fb = trial(b)
        if abs(fb) < best_res:
            best_res, best_dv = abs(fb), b

        grow = 0
        while fa * fb > 0 and grow < 80:
            b *= 1.6
            fb = trial(b)
            grow += 1
            if abs(fb) < best_res:
                best_res, best_dv = abs(fb), b
            if abs(b) > 50.0:  # unreasonable delta-v: give up the expansion
                break

        if fa * fb > 0:
            iterations = grow
        else:
            # Bisection: we rely only on the SIGN of the residual (not its
            # magnitude), which stays robust even when a bound falls into the
            # escape regime ("infinite" residual). Convergence if the residual
            # drops below the tolerance, OR if the delta-v interval becomes
            # negligible — in which case we hold the best possible solution
            # even for a very sensitive target (e.g. GEO apoapsis, where 1 m of
            # tolerance would represent a delta-v precision beyond numerical
            # reach).
            lo, flo = a, fa
            hi = b
            for iterations in range(1, max_iterations + 1):
                m = 0.5 * (lo + hi)
                fm = trial(m)
                if abs(fm) < best_res:
                    best_res, best_dv = abs(fm), m
                if abs(fm) < tolerance:
                    converged = True
                    best_dv, best_res = m, abs(fm)
                    break
                if abs(hi - lo) < 1e-10 * (1 + abs(m)):
                    converged = True
                    best_dv, best_res = m, abs(fm)
                    break
                if flo * fm < 0:
                    hi = m
                else:
                    lo, flo = m, fm

        # Return to the starting state before applying (or not) the real burn.
        self._restore_state(entry, snap)
        gator.UpdateSpaceObject()

        if not converged:
            return {
                "error": (
                    f"Targeting did not converge after {iterations} iterations "
                    f"(residual {best_res:.3e}). Best delta-v = {best_dv:.6f} km/s. "
                    f"Try a different initial_guess_km_s or check that the target "
                    f"is reachable with control '{control}'."
                ),
                "converged": False,
                "iterations": iterations,
                "best_delta_v_km_s": best_dv,
                "best_residual": best_res,
            }

        # Actual application of the converged burn (clean timeline + history).
        dv_components = [0.0, 0.0, 0.0]
        dv_components[axis - 1] = best_dv
        burn_result = self.apply_impulsive_burn(
            propagator_name, dv_components[0], dv_components[1],
            dv_components[2], frame="VNB",
        )
        achieved = measure(burn_result)

        return {
            "name": sat_name,
            "converged": True,
            "iterations": iterations,
            "goal_type": goal_type,
            "goal_value": goal_value,
            "achieved_value": achieved,
            "residual": achieved - goal_value,
            "delta_v_km_s": best_dv,
            "control": control,
            "burn_name": burn_result["burn_name"],
            "sma_km": burn_result["sma_km"],
            "ecc": burn_result["ecc"],
            "inc_deg": burn_result["inc_deg"],
            "raan_deg": burn_result["raan_deg"],
            "aop_deg": burn_result["aop_deg"],
            "ta_deg": burn_result["ta_deg"],
        }

    def time_to_apsis(self, propagator_name: str, apsis: str = "periapsis") -> dict:
        """
        Time remaining until the next periapsis/apoapsis/node from the
        propagator's current state (two-body Keplerian resolution).
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")
        valid = {"periapsis", "apoapsis", "ascending_node", "descending_node"}
        if apsis not in valid:
            raise ValueError(f"apsis '{apsis}' unknown — expected one of {sorted(valid)}.")

        entry = self.propagators[propagator_name]
        gator = entry["gator"]
        gator.UpdateSpaceObject()
        kep = _cartesian_to_keplerian(gator.GetState())
        sma = kep["sma_km"]
        ecc = kep["ecc"]
        if sma <= 0:
            raise ValueError("Unbound orbit (escape): time to apsis undefined.")

        mu = EARTH_MU_KM3_S2
        n = math.sqrt(mu / sma**3)  # mean motion (rad/s)
        period = 2 * math.pi / n
        two_pi = 2 * math.pi

        def ta_to_mean(ta_rad):
            E = 2 * math.atan2(math.sqrt(1 - ecc) * math.sin(ta_rad / 2),
                               math.sqrt(1 + ecc) * math.cos(ta_rad / 2))
            return (E - ecc * math.sin(E)) % two_pi

        ta_cur = math.radians(kep["ta_deg"])
        m_cur = ta_to_mean(ta_cur)

        if apsis == "periapsis":
            m_target = 0.0
        elif apsis == "apoapsis":
            m_target = math.pi
        else:
            # Nodes: argument of latitude u = aop + ta = 0 (ascending) or π.
            aop = math.radians(kep["aop_deg"])
            u_target = 0.0 if apsis == "ascending_node" else math.pi
            ta_target = (u_target - aop) % two_pi
            m_target = ta_to_mean(ta_target)

        dt = ((m_target - m_cur) % two_pi) / n

        return {
            "propagator_name": propagator_name,
            "apsis": apsis,
            "time_to_apsis_sec": dt,
            "period_sec": period,
            "current_ta_deg": kep["ta_deg"],
            "sma_km": sma,
            "ecc": ecc,
        }

    def get_history(self, propagator_name: str) -> list:
        """
        Returns the sampled state history (one point per Step() of
        propagate(), plus one post-burn point per maneuver), enriched with
        derived quantities: altitude, scalar speed, Keplerian elements.
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")

        entry = self.propagators[propagator_name]
        history = entry.get("history", [])
        if not history:
            raise ValueError(
                f"No history for '{propagator_name}' — call propagate() first."
            )

        samples = []
        for point in history:
            rx, ry, rz, vx, vy, vz = point["cartesian"]
            r = math.sqrt(rx**2 + ry**2 + rz**2)
            v = math.sqrt(vx**2 + vy**2 + vz**2)
            samples.append({
                "elapsed_sec": point["elapsed_sec"],
                "position_km": [rx, ry, rz],
                "velocity_km_s": [vx, vy, vz],
                "altitude_km": r - EARTH_RADIUS_KM,
                "speed_km_s": v,
                **_cartesian_to_keplerian(point["cartesian"]),
            })
        return samples

    def get_ground_track(self, propagator_name: str) -> list:
        """
        Ground track (geocentric latitude/longitude) computed from the state
        history, rotating the Earth at the GMST.

        Assumed approximation: uniform rotation at the sidereal rate from the
        GMST of the creation epoch (no nutation or polar motion, geocentric
        and not geodetic latitude) — consistent with the project's 2-body
        model, error ~ tenths of a degree.
        """
        if propagator_name not in self.propagators:
            raise ValueError(f"No propagator named '{propagator_name}' in this session.")

        entry = self.propagators[propagator_name]
        history = entry.get("history", [])
        if not history:
            raise ValueError(
                f"No history for '{propagator_name}' — call propagate() first."
            )

        sat_name = entry["spacecraft_name"]
        epoch = self.spacecraft_creation_state[sat_name]["epoch"]
        gmst0 = _gmst_rad(_epoch_to_jd(epoch))

        track = []
        for point in history:
            rx, ry, rz = point["cartesian"][:3]
            r = math.sqrt(rx**2 + ry**2 + rz**2)
            theta = gmst0 + EARTH_ROTATION_RAD_S * point["elapsed_sec"]
            lon = math.degrees(math.atan2(ry, rx) - theta)
            lon = (lon + 180.0) % 360.0 - 180.0  # normalize to ±180°
            lat = math.degrees(math.asin(rz / r))
            track.append({
                "elapsed_sec": point["elapsed_sec"],
                "lat_deg": lat,
                "lon_deg": lon,
            })
        return track

    def get_report(self) -> dict:
        satellites_report = []
        for name, sat in self.spacecraft.items():
            latest_state = None
            latest_source = None
            for entry in self.propagators.values():
                if entry["spacecraft_name"] == name and "last_keplerian_state" in entry:
                    latest_state = entry["last_keplerian_state"]
                    latest_source = entry.get("last_state_source", "propagated")

            if latest_state:
                satellites_report.append({
                    "name": name,
                    "epoch": sat.GetField("Epoch"),
                    "source": latest_source,
                    **latest_state,
                })
            else:
                satellites_report.append({
                    "name": name,
                    "epoch": sat.GetField("Epoch"),
                    "source": "creation",
                    "sma_km": float(sat.GetField("SMA")),
                    "ecc": float(sat.GetField("ECC")),
                    "inc_deg": float(sat.GetField("INC")),
                    "raan_deg": float(sat.GetField("RAAN")),
                    "aop_deg": float(sat.GetField("AOP")),
                    "ta_deg": float(sat.GetField("TA")),
                })

        propagators_report = []
        for prop_name, entry in self.propagators.items():
            maneuvers = [
                {"burn_name": t["burn_name"], "frame": t["frame"], "delta_v_km_s": t["delta_v_km_s"]}
                for t in entry.get("timeline", [])
                if t["type"] == "maneuver"
            ]
            propagators_report.append({
                "propagator_name": prop_name,
                "spacecraft_name": entry["spacecraft_name"],
                "step_size_sec": entry["step_size_sec"],
                "total_propagated_sec": entry.get("total_propagated_sec", 0.0),
                "maneuvers": maneuvers,
            })

        return {
            "satellites": satellites_report,
            "propagators": propagators_report,
            "satellite_count": len(satellites_report),
            "propagator_count": len(propagators_report),
        }

    def export_script(self, filename: str) -> dict:
        """
        Generates a GMAT .script file from the current session state
        (created spacecraft + propagators). Designed to be opened directly
        in desktop GMAT (GUI).

        Note: the API's SaveScript() was tested and does not work here
        -- it does not capture objects created via gmat.Construct(), so we
        rebuild the .script syntax ourselves.
        """
        if not self.spacecraft:
            raise ValueError("No spacecraft in the session, nothing to export.")

        lines = []
        lines.append("% GMAT script generated automatically by gmat-mcp")
        lines.append("% Generated via export_script() -- Phase 3")
        lines.append("")

        for name in self.spacecraft:
            # IMPORTANT: we declare the spacecraft in its CREATION state, not
            # its final state (last_keplerian_state) -- the Propagate/Maneuver
            # timeline below is what is meant to evolve the state from that
            # starting point. Using the final state here would replay the
            # entire mission on top of an already-arrived spacecraft (bug found
            # in practice while testing the export of a Hohmann transfer,
            # cf. Phase 5: the LEO departure orbit ended up declared as the GEO
            # arrival orbit).
            creation = self.spacecraft_creation_state[name]

            lines.append(f"Create Spacecraft {name};")
            lines.append(f"{name}.DateFormat = UTCGregorian;")
            lines.append(f"{name}.Epoch = '{creation['epoch']}';")
            lines.append(f"{name}.CoordinateSystem = EarthMJ2000Eq;")
            lines.append(f"{name}.DisplayStateType = Keplerian;")
            lines.append(f"{name}.SMA = {creation['sma_km']};")
            lines.append(f"{name}.ECC = {creation['ecc']};")
            lines.append(f"{name}.INC = {creation['inc_deg']};")
            lines.append(f"{name}.RAAN = {creation['raan_deg']};")
            lines.append(f"{name}.AOP = {creation['aop_deg']};")
            lines.append(f"{name}.TA = {creation['ta_deg']};")
            lines.append("")

        for prop_name, entry in self.propagators.items():
            sat_name = entry["spacecraft_name"]
            fm_name = f"FM_{sat_name}"
            fc = entry.get("force_config", {})

            # Spacecraft physical properties (declared here, before the
            # mission sequence — GMAT accepts them at this level).
            if fc.get("enable_drag") or fc.get("enable_srp"):
                lines.append(f"{sat_name}.DryMass = {fc['dry_mass_kg']};")
            if fc.get("enable_drag"):
                lines.append(f"{sat_name}.Cd = {fc['cd']};")
                lines.append(f"{sat_name}.DragArea = {fc['drag_area_m2']};")
            if fc.get("enable_srp"):
                lines.append(f"{sat_name}.Cr = {fc['cr']};")
                lines.append(f"{sat_name}.SRPArea = {fc['srp_area_m2']};")
            if fc.get("enable_drag") or fc.get("enable_srp"):
                lines.append("")

            lines.append(f"Create ForceModel {fm_name};")
            lines.append(f"{fm_name}.CentralBody = Earth;")

            # Gravity: harmonics or point-mass, + third-bodies.
            point_masses = ["Earth"] + list(fc.get("third_bodies", []))
            if fc.get("gravity_degree", 0) > 0:
                # The Earth is modeled by the harmonic field, not as a point-mass.
                point_masses = list(fc.get("third_bodies", []))
                if point_masses:
                    lines.append(f"{fm_name}.PointMasses = {{{', '.join(point_masses)}}};")
                lines.append(f"{fm_name}.GravityField.Earth.Degree = {fc['gravity_degree']};")
                lines.append(f"{fm_name}.GravityField.Earth.Order = {fc['gravity_order']};")
                lines.append(
                    f"{fm_name}.GravityField.Earth.PotentialFile = "
                    f"'{fc['gravity_file']}';"
                )
            else:
                lines.append(f"{fm_name}.PointMasses = {{{', '.join(point_masses)}}};")

            if fc.get("enable_drag"):
                lines.append(f"{fm_name}.Drag.AtmosphereModel = {fc['atmosphere_model']};")
                lines.append(f"{fm_name}.Drag.F107 = 150;")
                lines.append(f"{fm_name}.Drag.F107A = 150;")
                lines.append(f"{fm_name}.Drag.MagneticIndex = 3;")
            if fc.get("enable_srp"):
                lines.append(f"{fm_name}.SRP = On;")
            lines.append("")
            lines.append(f"Create Propagator {prop_name};")
            lines.append(f"{prop_name}.FM = {fm_name};")
            lines.append(f"{prop_name}.Type = PrinceDormand78;")
            lines.append(f"{prop_name}.InitialStepSize = {entry['step_size_sec']};")
            lines.append(f"{prop_name}.Accuracy = 1e-12;")
            lines.append(f"{prop_name}.MinStep = 0.0;")
            lines.append("")

            # One maneuver per distinct burn in the timeline (see Phase 5) --
            # same frame and components as those applied via apply_impulsive_burn().
            for t in entry.get("timeline", []):
                if t["type"] != "maneuver":
                    continue
                burn_name = t["burn_name"]
                dv1, dv2, dv3 = t["delta_v_km_s"]
                lines.append(f"Create ImpulsiveBurn {burn_name};")
                if t["frame"] == "Inertial":
                    lines.append(f"{burn_name}.CoordinateSystem = EarthMJ2000Eq;")
                else:
                    lines.append(f"{burn_name}.CoordinateSystem = Local;")
                    lines.append(f"{burn_name}.Origin = Earth;")
                    lines.append(f"{burn_name}.Axes = {t['frame']};")
                lines.append(f"{burn_name}.Element1 = {dv1};")
                lines.append(f"{burn_name}.Element2 = {dv2};")
                lines.append(f"{burn_name}.Element3 = {dv3};")
                lines.append("")

            # Finite-burn maneuvers: a dedicated tank + thruster + FiniteBurn
            # per burn (direction fixed on the thruster in the VNB frame).
            # Note: API-side dynamics discretized into impulses -> the result
            # replayed here in GMAT may differ by ~0.1 % (thrust arc).
            fb_list = [t for t in entry.get("timeline", []) if t["type"] == "finite_burn"]
            if fb_list:
                tank_names, thr_names = [], []
                for i, t in enumerate(fb_list, 1):
                    tk = f"{t['tank_name']}_{i}"
                    th = f"{t['thruster_name']}_{i}"
                    tank_names.append(tk)
                    thr_names.append(th)
                    init_fuel = self.thrusters.get(t["thruster_name"], {}).get(
                        "fuel_mass_kg", 0.0)
                    dx, dy, dz = t["direction_vec"]
                    lines.append(f"Create ChemicalTank {tk};")
                    lines.append(f"{tk}.FuelMass = {init_fuel};")
                    lines.append(f"Create ChemicalThruster {th};")
                    lines.append(f"{th}.C1 = {t['thrust_n']};")
                    lines.append(f"{th}.K1 = {t['isp_s']};")
                    lines.append(f"{th}.DecrementMass = true;")
                    lines.append(f"{th}.Tank = {{{tk}}};")
                    lines.append(f"{th}.CoordinateSystem = Local;")
                    lines.append(f"{th}.Origin = Earth;")
                    lines.append(f"{th}.Axes = VNB;")
                    lines.append(f"{th}.ThrustDirection1 = {dx};")
                    lines.append(f"{th}.ThrustDirection2 = {dy};")
                    lines.append(f"{th}.ThrustDirection3 = {dz};")
                    lines.append(f"Create FiniteBurn {t['burn_name']};")
                    lines.append(f"{t['burn_name']}.Thrusters = {{{th}}};")
                    lines.append("")
                lines.append(f"{sat_name}.Tanks = {{{', '.join(tank_names)}}};")
                lines.append(f"{sat_name}.Thrusters = {{{', '.join(thr_names)}}};")
                lines.append("")

        lines.append("BeginMissionSequence;")
        for prop_name, entry in self.propagators.items():
            sat_name = entry["spacecraft_name"]
            # Replays the timeline (propagations + maneuvers) in real
            # chronological order, rather than a single global Propagate --
            # necessary as soon as a maneuver has been inserted between two
            # propagations (cf. Phase 5).
            for t in entry.get("timeline", []):
                if t["type"] == "propagate":
                    lines.append(
                        f"Propagate {prop_name}({sat_name}) "
                        f"{{{sat_name}.ElapsedSecs = {t['duration_sec']}}};"
                    )
                elif t["type"] == "maneuver":
                    lines.append(f"Maneuver '{t['burn_name']}' {t['burn_name']}({sat_name});")
                elif t["type"] == "finite_burn":
                    bn = t["burn_name"]
                    lines.append(f"BeginFiniteBurn '{bn}' {bn}({sat_name});")
                    lines.append(
                        f"Propagate {prop_name}({sat_name}) "
                        f"{{{sat_name}.ElapsedSecs = {t['duration_sec']}}};"
                    )
                    lines.append(f"EndFiniteBurn '{bn}' {bn}({sat_name});")

        script_content = "\n".join(lines)
        filename = self.resolve_path(filename)

        with open(filename, "w") as f:
            f.write(script_content)

        return {
            "filename": filename,
            "satellites_exported": list(self.spacecraft.keys()),
            "propagators_exported": list(self.propagators.keys()),
        }


    def export_python(self, filename: str) -> dict:
        """
        Generates a standalone Python script that recreates the current
        session via the load_gmat/gmatpy API — reconstructs the creation
        state (not the propagation).
        """
        if not self.spacecraft:
            raise ValueError("No spacecraft in the session, nothing to export.")

        lines = []
        lines.append("# Python script generated automatically by gmat-mcp")
        lines.append("# Reproduces the session via the load_gmat/gmatpy API")
        lines.append("")
        lines.append("import load_gmat")
        lines.append("from load_gmat import gmat")
        lines.append("")

        for name, sat in self.spacecraft.items():
            var = name.lower()
            creation = self.spacecraft_creation_state[name]
            lines.append(f"{var} = gmat.Construct('Spacecraft', '{name}')")
            lines.append(f"{var}.SetField('DateFormat', 'UTCGregorian')")
            lines.append(f"{var}.SetField('Epoch', '{creation['epoch']}')")
            lines.append(f"{var}.SetField('CoordinateSystem', 'EarthMJ2000Eq')")
            lines.append(f"{var}.SetField('DisplayStateType', 'Keplerian')")
            lines.append(f"{var}.SetField('SMA', {creation['sma_km']})")
            lines.append(f"{var}.SetField('ECC', {creation['ecc']})")
            lines.append(f"{var}.SetField('INC', {creation['inc_deg']})")
            lines.append(f"{var}.SetField('RAAN', {creation['raan_deg']})")
            lines.append(f"{var}.SetField('AOP', {creation['aop_deg']})")
            lines.append(f"{var}.SetField('TA', {creation['ta_deg']})")

            # Physical properties if a propagator with drag/SRP uses them.
            for entry in self.propagators.values():
                if entry["spacecraft_name"] != name:
                    continue
                fc = entry.get("force_config", {})
                if fc.get("enable_drag") or fc.get("enable_srp"):
                    lines.append(f"{var}.SetField('DryMass', {fc['dry_mass_kg']})")
                if fc.get("enable_drag"):
                    lines.append(f"{var}.SetField('Cd', {fc['cd']})")
                    lines.append(f"{var}.SetField('DragArea', {fc['drag_area_m2']})")
                if fc.get("enable_srp"):
                    lines.append(f"{var}.SetField('Cr', {fc['cr']})")
                    lines.append(f"{var}.SetField('SRPArea', {fc['srp_area_m2']})")
                break
            lines.append("")

        # Reconstruction of the force models + propagators.
        for prop_name, entry in self.propagators.items():
            sat_name = entry["spacecraft_name"]
            fc = entry.get("force_config", {})
            fm_var = f"fm_{sat_name.lower()}"
            lines.append(f"{fm_var} = gmat.Construct('ForceModel', 'FM_{sat_name}')")
            lines.append(f"{fm_var}.SetField('CentralBody', 'Earth')")
            if fc.get("gravity_degree", 0) > 0:
                cof, _ = self._GRAVITY_FILES[fc["gravity_file"]]
                lines.append(f"_grav = gmat.Construct('GravityField')")
                lines.append(f"_grav.SetField('BodyName', 'Earth')")
                lines.append(f"_grav.SetField('Degree', {fc['gravity_degree']})")
                lines.append(f"_grav.SetField('Order', {fc['gravity_order']})")
                lines.append(f"_grav.SetField('PotentialFile', '{cof}')")
                lines.append(f"{fm_var}.AddForce(_grav)")
            else:
                lines.append(f"_pm = gmat.Construct('PointMassForce')")
                lines.append(f"_pm.SetField('BodyName', 'Earth')")
                lines.append(f"{fm_var}.AddForce(_pm)")
            for body in fc.get("third_bodies", []):
                lines.append(f"_tb = gmat.Construct('PointMassForce')")
                lines.append(f"_tb.SetField('BodyName', '{body}')")
                lines.append(f"{fm_var}.AddForce(_tb)")
            if fc.get("enable_drag"):
                lines.append(f"_drag = gmat.Construct('DragForce')")
                lines.append(f"_drag.SetField('AtmosphereModel', '{fc['atmosphere_model']}')")
                lines.append(f"_atmos = gmat.Construct('{fc['atmosphere_model']}')")
                lines.append(f"_drag.SetReference(_atmos)")
                lines.append(f"{fm_var}.AddForce(_drag)")
            if fc.get("enable_srp"):
                lines.append(f"_srp = gmat.Construct('SolarRadiationPressure')")
                lines.append(f"_srp.SetField('SRPModel', 'Spherical')")
                lines.append(f"{fm_var}.AddForce(_srp)")
            lines.append(
                f"{prop_name.lower()} = gmat.Construct('Propagator', '{prop_name}')"
            )
            lines.append(
                f"_integ = gmat.Construct('PrinceDormand78', 'Integ_{sat_name}')"
            )
            lines.append(f"{prop_name.lower()}.SetReference(_integ)")
            lines.append(f"{prop_name.lower()}.SetReference({fm_var})")
            lines.append(
                f"{prop_name.lower()}.SetField('InitialStepSize', {entry['step_size_sec']})"
            )
            lines.append(f"{prop_name.lower()}.SetField('Accuracy', 1e-12)")
            lines.append("")

        lines.append("gmat.Initialize()")
        lines.append("")
        lines.append("print('Session reconstructed successfully')")

        script_content = "\n".join(lines)
        filename = self.resolve_path(filename)

        with open(filename, "w") as f:
            f.write(script_content)

        return {
            "filename": filename,
            "satellites_exported": list(self.spacecraft.keys()),
        }
# Single instance shared across the whole server (singleton).
session = GmatSession()