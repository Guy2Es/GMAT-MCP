# server.py
from mcp.server.fastmcp import FastMCP
from tools.spacecraft import register_spacecraft_tools
from tools.propagation import register_propagation_tools
from tools.reports import register_report_tools
from tools.mission_planner import register_mission_planner_tools
from tools.maneuvers import register_maneuver_tools
from tools.analysis import register_analysis_tools
from tools.session import register_session_tools
from tools.targeting import register_targeting_tools
from tools.tle import register_tle_tools
from tools.astro_utils import register_astro_tools

mcp = FastMCP("gmat-mcp")

register_session_tools(mcp)
register_spacecraft_tools(mcp)
register_propagation_tools(mcp)
register_report_tools(mcp)
register_mission_planner_tools(mcp)
register_maneuver_tools(mcp)
register_analysis_tools(mcp)
register_targeting_tools(mcp)
register_tle_tools(mcp)
register_astro_tools(mcp)

if __name__ == "__main__":
    mcp.run()