# tools/targeting.py
"""
MCP targeting tools (in-house differential corrector): find the delta-v that
drives an orbital quantity to a target value.
"""

from gmat_session import session


def register_targeting_tools(mcp):

    @mcp.tool()
    def target_maneuver(
        propagator_name: str,
        goal_type: str,
        goal_value: float,
        control: str = "prograde",
        initial_guess_km_s: float = 0.1,
        tolerance: float = 1e-3,
        max_iterations: int = 60,
    ) -> dict:
        """
        Automatically searches for the delta-v (by iterative shooting / secant
        method) that drives a final orbital quantity to a target value, then
        applies the converged maneuver. Equivalent to a GMAT
        Target/Vary/Achieve sequence, in the VNB frame.

        Example: from a circular LEO orbit, target
        goal_type="apoapsis_radius_km", goal_value=42164 to find the prograde
        delta-v of the first stage of a Hohmann transfer.

        Args:
            propagator_name: Name of the propagator (created via create_propagator).
            goal_type: Quantity to target — "apoapsis_radius_km",
                "periapsis_radius_km", "apoapsis_altitude_km",
                "periapsis_altitude_km", "sma_km", "ecc" or "period_sec".
            goal_value: Target value (km, dimensionless, or seconds).
            control: Delta-v control axis — "prograde" (default, raises the
                apoapsis), "normal" or "binormal" (VNB frame).
            initial_guess_km_s: First delta-v trial (km/s).
            tolerance: Tolerance on the residual (same unit as goal_value).
            max_iterations: Maximum number of iterations (default 60).
        """
        try:
            return session.target_maneuver(
                propagator_name, goal_type, goal_value, control=control,
                initial_guess_km_s=initial_guess_km_s, tolerance=tolerance,
                max_iterations=max_iterations,
            )
        except Exception as e:
            return {"error": str(e)}
