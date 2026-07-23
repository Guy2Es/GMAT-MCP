# tools/session.py
"""
MCP session-management tools: mission name and output folder.
"""

from gmat_session import session


def register_session_tools(mcp):

    @mcp.tool()
    def set_mission(name: str, base_path: str = None) -> dict:
        """
        Sets the current mission name. All subsequent exports
        (export_script, export_python, export_data, plot_*) will
        automatically go into <base_path or project folder>/output/<name>/,
        without having to repeat the mission name on every call.

        Call once at the start of the session, with a name derived from the
        user request (e.g. "Hohmann_LEO_GEO", "Transfert_Terre_Lune").

        Args:
            name: Mission name (used as folder name).
            base_path: Optional base folder (default: project folder).
        """
        try:
            return session.set_mission(name, base_path)
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def get_mission_info() -> dict:
        """Returns the currently active mission name and output folder (if set)."""
        return {"mission_name": session.mission_name, "output_dir": session.output_dir}

    @mcp.tool()
    def reset_session() -> dict:
        """
        Resets the GMAT session: removes all spacecraft, propagators, force
        models and maneuvers to start from scratch, without restarting the
        server. The mission name and output folder are kept.

        Useful to reuse an already-taken name, or to chain a new mission
        within the same process.
        """
        try:
            return session.reset()
        except Exception as e:
            return {"error": str(e)}
