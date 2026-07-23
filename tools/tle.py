# tools/tle.py
"""
MCP tools to import TLEs (Two-Line Elements) via SGP4.
"""

from gmat_session import session


def register_tle_tools(mcp):

    @mcp.tool()
    def parse_tle(line1: str, line2: str) -> dict:
        """
        Reads a TLE (2 lines) and returns its mean elements (NORAD number,
        epoch, inclination, RAAN, eccentricity, argument of periapsis, mean
        anomaly, mean motion, derived semi-major axis, period, B*), without
        creating anything in GMAT.

        Args:
            line1: First line of the TLE (starts with "1 ").
            line2: Second line of the TLE (starts with "2 ").
        """
        try:
            return session.parse_tle(line1, line2)
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def create_spacecraft_from_tle(name: str, line1: str, line2: str) -> dict:
        """
        Creates a spacecraft from a TLE: propagates SGP4 to the TLE epoch,
        converts the TEME state to EarthMJ2000Eq (via GMAT) and instantiates
        the spacecraft at that inertial cartesian state.

        Call BEFORE create_propagator. Note: the TLE provides MEAN elements;
        the re-read osculating SMA may differ by ~10 km, which is expected.

        Args:
            name: Unique spacecraft name.
            line1: First line of the TLE.
            line2: Second line of the TLE.
        """
        try:
            return session.create_spacecraft_from_tle(name, line1, line2)
        except Exception as e:
            return {"error": str(e)}
