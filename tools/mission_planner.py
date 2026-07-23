# tools/mission_planner.py
"""
MCP orbital-calculation tools (two-body) complementing the space-mechanics
skill — no call to the GMAT API or the session, just formulas, to eliminate
the risk of arithmetic error on the most common calculations (Hohmann,
escape velocity, feasibility).
"""
import math

from gmat_session import EARTH_MU_KM3_S2, EARTH_RADIUS_KM


def register_mission_planner_tools(mcp):

    @mcp.tool()
    def suggest_hohmann_transfer(from_sma_km: float, to_sma_km: float) -> dict:
        """
        Computes a Hohmann transfer between two coplanar circular orbits
        (radii from the center of the Earth, not altitudes).

        Args:
            from_sma_km: Radius of the departure circular orbit, in km.
            to_sma_km: Radius of the arrival circular orbit, in km.
        """
        r1, r2 = from_sma_km, to_sma_km
        mu = EARTH_MU_KM3_S2

        a_transfer = (r1 + r2) / 2

        v_circ1 = math.sqrt(mu / r1)
        v_transfer_p = math.sqrt(mu * (2 / r1 - 1 / a_transfer))
        delta_v1 = abs(v_transfer_p - v_circ1)

        v_circ2 = math.sqrt(mu / r2)
        v_transfer_a = math.sqrt(mu * (2 / r2 - 1 / a_transfer))
        delta_v2 = abs(v_circ2 - v_transfer_a)

        transfer_time_sec = math.pi * math.sqrt(a_transfer**3 / mu)

        return {
            "from_sma_km": r1,
            "to_sma_km": r2,
            "transfer_sma_km": a_transfer,
            "delta_v1_km_s": delta_v1,
            "delta_v2_km_s": delta_v2,
            "delta_v_total_km_s": delta_v1 + delta_v2,
            "transfer_time_sec": transfer_time_sec,
        }

    @mcp.tool()
    def calculate_escape_velocity(alt_km: float) -> dict:
        """
        Computes the escape velocity and the reference circular velocity at
        a given altitude above the Earth.

        Args:
            alt_km: Altitude above the Earth's surface, in km.
        """
        r = EARTH_RADIUS_KM + alt_km
        mu = EARTH_MU_KM3_S2

        return {
            "alt_km": alt_km,
            "radius_km": r,
            "escape_velocity_km_s": math.sqrt(2 * mu / r),
            "circular_velocity_km_s": math.sqrt(mu / r),
        }

    @mcp.tool()
    def validate_orbit(sma_km: float, ecc: float, inc_deg: float) -> dict:
        """
        Checks the simple physical feasibility of an orbit (basic Keplerian
        guardrails: periapsis above the atmosphere, closed-orbit
        eccentricity, valid inclination).

        Args:
            sma_km: Semi-major axis, in km.
            ecc: Orbit eccentricity.
            inc_deg: Inclination, in degrees.
        """
        issues = []
        warnings = []

        perigee_km = sma_km * (1 - ecc)
        perigee_alt_km = perigee_km - EARTH_RADIUS_KM

        if perigee_alt_km <= 0:
            issues.append(
                f"Periapsis ({perigee_alt_km:.1f} km altitude) below the Earth's surface — impossible orbit."
            )
        elif perigee_alt_km < 150:
            warnings.append(
                f"Low periapsis ({perigee_alt_km:.1f} km altitude) — the orbit decays quickly (drag)."
            )

        if not (0 <= ecc < 1):
            issues.append(f"ECC={ecc} out of [0, 1) — non-closed orbit, not handled by the project propagator.")

        if not (0 <= inc_deg <= 180):
            issues.append(f"INC={inc_deg}° out of [0, 180°].")

        return {
            "sma_km": sma_km,
            "ecc": ecc,
            "inc_deg": inc_deg,
            "perigee_alt_km": perigee_alt_km,
            "feasible": len(issues) == 0,
            "issues": issues,
            "warnings": warnings,
        }
