# tools/reports.py
"""
MCP tools related to session reports and exports.
"""

from gmat_session import session


def register_report_tools(mcp):

    @mcp.tool()
    def get_report() -> dict:
        """
        Returns a structured summary of everything created in the current
        session: spacecraft (with their current orbital elements, propagated
        if applicable), and propagators.
        """
        try:
            return session.get_report()
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def export_script(filename: str) -> dict:
        """
        Exports the current session state as a GMAT .script file, directly
        openable in desktop GMAT (full GUI).

        If set_mission() was called, only the file name is used (no path) —
        the file goes into the current mission's folder.

        Args:
            filename: Name of the .script file to create (e.g. "mission.script").
        """
        try:
            return session.export_script(filename)
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def export_python(filename: str) -> dict:
        """
        Exports the current session state as a standalone Python script,
        reproducing the creation of the spacecraft via the load_gmat API.

        If set_mission() was called, only the file name is used (no path) —
        the file goes into the current mission's folder.

        Args:
            filename: Name of the .py file to create (e.g. "mission.py").
        """
        try:
            return session.export_python(filename)
        except Exception as e:
            return {"error": str(e)}