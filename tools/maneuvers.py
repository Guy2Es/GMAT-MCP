# tools/maneuvers.py
"""
MCP tools related to impulsive orbital maneuvers.
"""

from gmat_session import session


def register_maneuver_tools(mcp):

    @mcp.tool()
    def apply_impulsive_burn(
        propagator_name: str,
        delta_v1_km_s: float,
        delta_v2_km_s: float = 0.0,
        delta_v3_km_s: float = 0.0,
        frame: str = "VNB",
    ) -> dict:
        """
        Applies an impulsive maneuver (instantaneous delta-v, at unchanged
        position) to the spacecraft associated with this propagator, and
        returns the resulting orbital state.

        Args:
            propagator_name: Name of the propagator (created via create_propagator).
            delta_v1_km_s: First component of the delta-v, in km/s.
            delta_v2_km_s: Second component of the delta-v, in km/s.
            delta_v3_km_s: Third component of the delta-v, in km/s.
            frame: Delta-v frame -- "VNB" (prograde/normal/binormal,
                default; positive Element1 = prograde), "LVLH", "MJ2000Eq"
                (these three in an Earth-centered local frame), or "Inertial"
                (delta-v expressed directly in EarthMJ2000Eq).
        """
        try:
            return session.apply_impulsive_burn(
                propagator_name, delta_v1_km_s, delta_v2_km_s, delta_v3_km_s, frame
            )
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def create_thruster(
        spacecraft_name: str,
        thruster_name: str = None,
        thrust_n: float = 440.0,
        isp_s: float = 320.0,
        fuel_mass_kg: float = 200.0,
        tank_name: str = None,
    ) -> dict:
        """
        Declares a chemical thruster (constant thrust) and its tank for a
        spacecraft, as a prerequisite to a finite-burn maneuver
        (apply_finite_burn).

        Args:
            spacecraft_name: Name of the spacecraft that will carry the thruster.
            thruster_name: Thruster name (default "Thruster_<sat>").
            thrust_n: Thrust in newtons (default 440 N, typical apogee motor).
            isp_s: Specific impulse in seconds (default 320 s).
            fuel_mass_kg: Available propellant mass in kg (default 200 kg).
            tank_name: Tank name (default "Tank_<sat>").
        """
        try:
            return session.create_thruster(
                spacecraft_name, thruster_name=thruster_name,
                thrust_n=thrust_n, isp_s=isp_s,
                fuel_mass_kg=fuel_mass_kg, tank_name=tank_name,
            )
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def apply_finite_burn(
        propagator_name: str,
        duration_sec: float,
        direction: str = "prograde",
        thruster_name: str = None,
        segment_sec: float = None,
    ) -> dict:
        """
        Applies a finite-burn maneuver (continuous thrust) over a given
        duration, integrating the thrust segment by segment. Requires a
        thruster created beforehand via create_thruster.

        Returns the final orbital state, the total accumulated delta-v, the
        propellant consumed/remaining and the final mass.

        Args:
            propagator_name: Name of the propagator (created via create_propagator).
            duration_sec: Burn duration in seconds.
            direction: "prograde" (default), "retrograde", "normal" or
                "antinormal" (local VNB frame).
            thruster_name: Thruster to use (default: the spacecraft's only
                thruster).
            segment_sec: Duration of a discretization segment (default:
                ~ duration/20, capped by the propagator step).
        """
        try:
            return session.apply_finite_burn(
                propagator_name, duration_sec, direction=direction,
                thruster_name=thruster_name, segment_sec=segment_sec,
            )
        except Exception as e:
            return {"error": str(e)}
