# tools/spacecraft.py
"""
MCP tools related to creating and inspecting GMAT spacecraft.
"""

from gmat_session import session


def register_spacecraft_tools(mcp):
    """
    Registers the spacecraft tools on the FastMCP instance passed
    as argument. Called from server.py.
    """

    @mcp.tool()
    def create_spacecraft(
        name: str,
        sma_km: float = None,
        ecc: float = None,
        inc_deg: float = None,
        raan_deg: float = 0.0,
        aop_deg: float = 0.0,
        ta_deg: float = 0.0,
        epoch: str = "01 Jan 2026 00:00:00.000",
        cartesian_state: list = None,
    ) -> dict:
        """
        Creates a spacecraft in GMAT, either from Keplerian orbital elements
        or from an inertial cartesian state. The two modes are exclusive.

        Args:
            name: Unique spacecraft name (e.g. "Sat1").
            sma_km: Semi-major axis in kilometers (e.g. 6978 for ~600km altitude).
            ecc: Orbit eccentricity (0 = circular, between 0 and 1).
            inc_deg: Orbital inclination in degrees (0 = equatorial, 90 = polar).
            raan_deg: Right ascension of the ascending node, in degrees.
            aop_deg: Argument of periapsis, in degrees.
            ta_deg: Initial true anomaly, in degrees.
            epoch: UTCGregorian epoch (e.g. "01 Jan 2026 00:00:00.000").
            cartesian_state: Alternative to Keplerian elements — inertial state
                [x, y, z, vx, vy, vz] in km and km/s (EarthMJ2000Eq frame).
        """
        try:
            return session.create_spacecraft(
                name, sma_km, ecc, inc_deg, raan_deg, aop_deg, ta_deg,
                epoch=epoch, cartesian_state=cartesian_state,
            )
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def get_spacecraft_state(name: str) -> dict:
        """
        Re-reads the current orbital state of an already-created spacecraft.

        Args:
            name: Name of the spacecraft to query.
        """
        try:
            return session.get_spacecraft_state(name)
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def list_spacecraft() -> dict:
        """Lists all spacecraft created in the current GMAT session."""
        return {"spacecraft": session.list_spacecraft()}