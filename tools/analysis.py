# tools/analysis.py
"""
MCP visualization and data-export tools (Phase 6): 3D trajectories,
orbital-element evolution, ground track, CSV/JSON export.

The matplotlib plots are returned as a FastMCP Image (inline display in the
Claude Desktop chat) AND saved as PNG on disk; the interactive_html=True
parameter additionally produces an interactive plotly .html.
"""
import csv
import json
import math
import os

import matplotlib
matplotlib.use("Agg")  # headless server process -- must be done before pyplot
import matplotlib.pyplot as plt

from mcp.server.fastmcp import Image

from gmat_session import session, EARTH_RADIUS_KM


def _maneuver_points(propagator_name: str) -> list:
    """Returns [(burn_name, sample), ...] to mark the burns on the plots."""
    entry = session.propagators[propagator_name]
    history = entry.get("history", [])
    points = []
    for t in entry.get("timeline", []):
        if t["type"] != "maneuver" or "at_elapsed_sec" not in t:
            continue
        # The post-burn sample is the last one in the history at that instant.
        sample = None
        for point in history:
            if point["elapsed_sec"] <= t["at_elapsed_sec"]:
                sample = point
        if sample:
            points.append((t["burn_name"], sample))
    return points


def _time_axis(samples: list):
    """Readable time axis: seconds below 2 h, hours beyond."""
    max_elapsed = samples[-1]["elapsed_sec"]
    if max_elapsed < 7200:
        return [s["elapsed_sec"] for s in samples], "Elapsed time (s)"
    return [s["elapsed_sec"] / 3600 for s in samples], "Elapsed time (h)"


def register_analysis_tools(mcp):

    @mcp.tool()
    def export_data(propagator_name: str, filename: str, format: str = "csv") -> dict:
        """
        Exports the propagated state history (time, position, velocity,
        altitude, Keplerian elements) as CSV or JSON for external analysis
        (Excel, matplotlib, pandas...).

        Args:
            propagator_name: Name of the propagator (already propagated at least once).
            filename: Output file (.csv or .json). If set_mission() was called,
                only the file name is used (no path) — the file goes into the
                current mission's folder; otherwise, relative = project folder.
            format: "csv" (default) or "json".
        """
        try:
            if format not in ("csv", "json"):
                raise ValueError(f"Format '{format}' unknown — expected 'csv' or 'json'.")

            samples = session.get_history(propagator_name)
            path = session.resolve_path(filename)

            if format == "json":
                with open(path, "w") as f:
                    json.dump(samples, f, indent=2)
            else:
                columns = ["elapsed_sec",
                           "rx_km", "ry_km", "rz_km",
                           "vx_km_s", "vy_km_s", "vz_km_s",
                           "altitude_km", "speed_km_s",
                           "sma_km", "ecc", "inc_deg",
                           "raan_deg", "aop_deg", "ta_deg"]
                with open(path, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(columns)
                    for s in samples:
                        writer.writerow([
                            s["elapsed_sec"],
                            *s["position_km"], *s["velocity_km_s"],
                            s["altitude_km"], s["speed_km_s"],
                            s["sma_km"], s["ecc"], s["inc_deg"],
                            s["raan_deg"], s["aop_deg"], s["ta_deg"],
                        ])

            return {
                "filename": path,
                "format": format,
                "sample_count": len(samples),
                "duration_sec": samples[-1]["elapsed_sec"] - samples[0]["elapsed_sec"],
            }
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def plot_trajectory_3d(propagator_names: list = None,
                           filename: str = "trajectory_3d.png",
                           interactive_html: bool = False):
        """
        Plots the orbital trajectory (or trajectories) in 3D around the
        Earth, with markers at the maneuver points. The image is shown in the
        chat and the PNG is saved on disk.

        Args:
            propagator_names: List of propagators to plot (default: all those
                that have already been propagated).
            filename: Output PNG file. If set_mission() was called, only the
                file name is used — the file goes into the current mission's
                folder; otherwise, relative = project folder.
            interactive_html: If True, additionally writes an interactive
                plotly .html (mouse rotation) next to the PNG.
        """
        try:
            if propagator_names is None:
                propagator_names = [name for name, entry in session.propagators.items()
                                    if entry.get("history")]
            if not propagator_names:
                raise ValueError("No propagator with a history — call propagate() first.")

            all_series = {name: session.get_history(name) for name in propagator_names}

            fig = plt.figure(figsize=(9, 9))
            ax = fig.add_subplot(projection="3d")

            # Earth sphere at the correct radius
            import numpy as np
            u = np.linspace(0, 2 * np.pi, 40)
            v = np.linspace(0, np.pi, 20)
            xs = EARTH_RADIUS_KM * np.outer(np.cos(u), np.sin(v))
            ys = EARTH_RADIUS_KM * np.outer(np.sin(u), np.sin(v))
            zs = EARTH_RADIUS_KM * np.outer(np.ones_like(u), np.cos(v))
            ax.plot_surface(xs, ys, zs, color="steelblue", alpha=0.4, linewidth=0)

            max_range = EARTH_RADIUS_KM
            for name, samples in all_series.items():
                sat_name = session.propagators[name]["spacecraft_name"]
                x = [s["position_km"][0] for s in samples]
                y = [s["position_km"][1] for s in samples]
                z = [s["position_km"][2] for s in samples]
                ax.plot(x, y, z, label=sat_name, linewidth=1.2)
                max_range = max(max_range, *(abs(c) for c in x + y + z))

                for burn_name, sample in _maneuver_points(name):
                    px, py, pz = sample["cartesian"][:3]
                    ax.scatter([px], [py], [pz], color="red", s=45, marker="^", zorder=5)
                    ax.text(px, py, pz, f"  {burn_name}", fontsize=7, color="darkred")

            # Equal-scale axes (otherwise the orbit is distorted)
            ax.set_xlim(-max_range, max_range)
            ax.set_ylim(-max_range, max_range)
            ax.set_zlim(-max_range, max_range)
            ax.set_box_aspect([1, 1, 1])
            ax.set_xlabel("X (km)")
            ax.set_ylabel("Y (km)")
            ax.set_zlabel("Z (km)")
            ax.set_title("Orbital trajectories (EarthMJ2000Eq)")
            ax.legend(loc="upper right", fontsize=8)

            png_path = session.resolve_path(filename)
            fig.savefig(png_path, dpi=110, bbox_inches="tight")
            plt.close(fig)

            html_path = None
            if interactive_html:
                import plotly.graph_objects as go
                pfig = go.Figure()
                pfig.add_surface(x=xs, y=ys, z=zs, opacity=0.4, showscale=False,
                                 colorscale=[[0, "steelblue"], [1, "steelblue"]])
                for name, samples in all_series.items():
                    sat_name = session.propagators[name]["spacecraft_name"]
                    pfig.add_scatter3d(
                        x=[s["position_km"][0] for s in samples],
                        y=[s["position_km"][1] for s in samples],
                        z=[s["position_km"][2] for s in samples],
                        mode="lines", name=sat_name)
                    burns = _maneuver_points(name)
                    if burns:
                        pfig.add_scatter3d(
                            x=[s["cartesian"][0] for _, s in burns],
                            y=[s["cartesian"][1] for _, s in burns],
                            z=[s["cartesian"][2] for _, s in burns],
                            mode="markers+text",
                            text=[b for b, _ in burns],
                            marker=dict(color="red", size=5, symbol="diamond"),
                            name="Maneuvers")
                pfig.update_layout(
                    title="Orbital trajectories (EarthMJ2000Eq)",
                    scene=dict(aspectmode="data",
                               xaxis_title="X (km)", yaxis_title="Y (km)",
                               zaxis_title="Z (km)"))
                html_path = os.path.splitext(png_path)[0] + ".html"
                pfig.write_html(html_path, include_plotlyjs=True)

            summary = {"png": png_path, "html": html_path,
                       "propagators_plotted": propagator_names}
            return [Image(path=png_path), json.dumps(summary)]
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def plot_orbital_elements(propagator_name: str,
                              filename: str = None,
                              interactive_html: bool = False):
        """
        Plots the time evolution of the orbital elements (SMA, ECC, INC,
        altitude) with vertical lines at the maneuvers. The image is shown in
        the chat and the PNG is saved on disk.

        Args:
            propagator_name: Name of the propagator (already propagated at least once).
            filename: Output PNG file (default: orbital_elements_<name>.png).
                If set_mission() was called, the file goes into the current
                mission's folder.
            interactive_html: If True, additionally writes an interactive plotly .html.
        """
        try:
            samples = session.get_history(propagator_name)
            t, t_label = _time_axis(samples)
            t_scale = 3600 if t_label.endswith("(h)") else 1
            burns = [(b, s["elapsed_sec"] / t_scale) for b, s in _maneuver_points(propagator_name)]

            series = [
                ("SMA (km)", [s["sma_km"] for s in samples]),
                ("ECC", [s["ecc"] for s in samples]),
                ("INC (deg)", [s["inc_deg"] for s in samples]),
                ("Altitude (km)", [s["altitude_km"] for s in samples]),
            ]

            fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
            sat_name = session.propagators[propagator_name]["spacecraft_name"]
            fig.suptitle(f"Orbital elements — {sat_name}")
            for ax, (label, values) in zip(axes.flat, series):
                ax.plot(t, values, linewidth=1.2)
                for burn_name, burn_t in burns:
                    ax.axvline(burn_t, color="red", linestyle="--", linewidth=0.8, alpha=0.7)
                ax.set_ylabel(label)
                ax.grid(True, alpha=0.3)
            for ax in axes[-1]:
                ax.set_xlabel(t_label)

            png_path = session.resolve_path(filename or f"orbital_elements_{propagator_name}.png")
            fig.savefig(png_path, dpi=110, bbox_inches="tight")
            plt.close(fig)

            html_path = None
            if interactive_html:
                import plotly.graph_objects as go
                from plotly.subplots import make_subplots
                pfig = make_subplots(rows=2, cols=2,
                                     subplot_titles=[label for label, _ in series])
                for i, (label, values) in enumerate(series):
                    pfig.add_scatter(x=t, y=values, mode="lines", name=label,
                                     row=i // 2 + 1, col=i % 2 + 1)
                for burn_name, burn_t in burns:
                    pfig.add_vline(x=burn_t, line_dash="dash", line_color="red")
                pfig.update_layout(title=f"Orbital elements — {sat_name}",
                                   showlegend=False)
                pfig.update_xaxes(title_text=t_label, row=2)
                html_path = os.path.splitext(png_path)[0] + ".html"
                pfig.write_html(html_path, include_plotlyjs=True)

            summary = {"png": png_path, "html": html_path,
                       "maneuvers_marked": [b for b, _ in burns]}
            return [Image(path=png_path), json.dumps(summary)]
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def plot_ground_track(propagator_name: str,
                          filename: str = None,
                          interactive_html: bool = False):
        """
        Plots the spacecraft's ground track (latitude/longitude), computed
        with the Earth's real rotation (GMST). The image is shown in the chat
        and the PNG is saved on disk.

        Args:
            propagator_name: Name of the propagator (already propagated at least once).
            filename: Output PNG file (default: ground_track_<name>.png).
                If set_mission() was called, the file goes into the current
                mission's folder.
            interactive_html: If True, additionally writes an interactive plotly .html.
        """
        try:
            track = session.get_ground_track(propagator_name)
            sat_name = session.propagators[propagator_name]["spacecraft_name"]

            # Segments cut at the date-line crossing, otherwise matplotlib
            # draws a spurious horizontal line at ±180°.
            segments = [[]]
            for i, p in enumerate(track):
                if i > 0 and abs(p["lon_deg"] - track[i - 1]["lon_deg"]) > 180:
                    segments.append([])
                segments[-1].append(p)

            fig, ax = plt.subplots(figsize=(11, 5.5))
            for seg in segments:
                ax.plot([p["lon_deg"] for p in seg], [p["lat_deg"] for p in seg],
                        color="tab:blue", linewidth=1.0)
            ax.plot(track[0]["lon_deg"], track[0]["lat_deg"], "go", markersize=7,
                    label="Start")
            ax.plot(track[-1]["lon_deg"], track[-1]["lat_deg"], "rs", markersize=7,
                    label="End")
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            ax.set_xticks(range(-180, 181, 30))
            ax.set_yticks(range(-90, 91, 30))
            ax.grid(True, alpha=0.4)
            ax.set_xlabel("Longitude (deg)")
            ax.set_ylabel("Latitude (deg)")
            ax.set_title(f"Ground track — {sat_name}")
            ax.legend(loc="upper right", fontsize=8)

            png_path = session.resolve_path(filename or f"ground_track_{propagator_name}.png")
            fig.savefig(png_path, dpi=110, bbox_inches="tight")
            plt.close(fig)

            html_path = None
            if interactive_html:
                import plotly.graph_objects as go
                pfig = go.Figure()
                lons, lats = [], []
                for seg in segments:
                    lons += [p["lon_deg"] for p in seg] + [None]
                    lats += [p["lat_deg"] for p in seg] + [None]
                pfig.add_scatter(x=lons, y=lats, mode="lines", name=sat_name)
                pfig.update_layout(
                    title=f"Ground track — {sat_name}",
                    xaxis=dict(title="Longitude (deg)", range=[-180, 180], dtick=30),
                    yaxis=dict(title="Latitude (deg)", range=[-90, 90], dtick=30))
                html_path = os.path.splitext(png_path)[0] + ".html"
                pfig.write_html(html_path, include_plotlyjs=True)

            summary = {"png": png_path, "html": html_path,
                       "point_count": len(track)}
            return [Image(path=png_path), json.dumps(summary)]
        except Exception as e:
            return {"error": str(e)}
