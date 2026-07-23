# tools/astro_utils.py
"""
Analytic MCP tools (pure Python, no GMAT): transfers, plane change,
Lambert, maneuver timing, time conversions, equinoctial elements,
Tsiolkovsky, gravitational parameters.
"""
import math

from gmat_session import session, EARTH_MU_KM3_S2, EARTH_RADIUS_KM

# Standard gravitational parameters (km^3/s^2) and equatorial radii (km).
MU = {
    "Sun": 1.32712440018e11,
    "Mercury": 22032.09,
    "Venus": 324858.592,
    "Earth": EARTH_MU_KM3_S2,
    "Moon": 4902.8005821478,
    "Mars": 42828.375214,
    "Jupiter": 126712764.8,
    "Saturn": 37940585.2,
    "Uranus": 5794548.6,
    "Neptune": 6836527.10058,
}
BODY_RADIUS_KM = {
    "Sun": 696000.0, "Mercury": 2439.7, "Venus": 6051.8, "Earth": EARTH_RADIUS_KM,
    "Moon": 1737.4, "Mars": 3396.19, "Jupiter": 71492.0, "Saturn": 60268.0,
    "Uranus": 25559.0, "Neptune": 24764.0,
}
G0 = 9.80665  # m/s^2, for Tsiolkovsky


def _stumpff(psi):
    if psi > 1e-6:
        sq = math.sqrt(psi)
        c2 = (1 - math.cos(sq)) / psi
        c3 = (sq - math.sin(sq)) / (sq**3)
    elif psi < -1e-6:
        sq = math.sqrt(-psi)
        c2 = (1 - math.cosh(sq)) / psi
        c3 = (math.sinh(sq) - sq) / (sq**3)
    else:
        c2, c3 = 0.5, 1.0 / 6.0
    return c2, c3


def register_astro_tools(mcp):

    @mcp.tool()
    def suggest_bielliptic_transfer(r1_km: float, r2_km: float, rb_km: float) -> dict:
        """
        Computes a bi-elliptic transfer (3 impulses) between two coplanar
        circular orbits of radii r1 and r2, passing through an intermediate
        radius rb (common apoapsis, rb >= r2). Automatically compares it to
        the Hohmann transfer to tell which is more economical.

        Args:
            r1_km: Radius of the departure orbit (km, from the center).
            r2_km: Radius of the arrival orbit (km).
            rb_km: Intermediate apoapsis radius (km, typically > r2).
        """
        try:
            mu = EARTH_MU_KM3_S2
            if rb_km < r2_km:
                return {"error": "rb_km must be >= r2_km for a bi-elliptic transfer."}
            a1 = (r1_km + rb_km) / 2.0
            a2 = (rb_km + r2_km) / 2.0
            vc1 = math.sqrt(mu / r1_km)
            vc2 = math.sqrt(mu / r2_km)
            vp1 = math.sqrt(mu * (2.0 / r1_km - 1.0 / a1))
            va1 = math.sqrt(mu * (2.0 / rb_km - 1.0 / a1))
            vp2 = math.sqrt(mu * (2.0 / rb_km - 1.0 / a2))
            va2 = math.sqrt(mu * (2.0 / r2_km - 1.0 / a2))
            dv1 = vp1 - vc1
            dv2 = vp2 - va1
            dv3 = va2 - vc2  # negative = braking at arrival
            total = abs(dv1) + abs(dv2) + abs(dv3)

            # Reference Hohmann for comparison.
            ah = (r1_km + r2_km) / 2.0
            vph = math.sqrt(mu * (2.0 / r1_km - 1.0 / ah))
            vah = math.sqrt(mu * (2.0 / r2_km - 1.0 / ah))
            hohmann_total = abs(vph - vc1) + abs(vc2 - vah)

            return {
                "delta_v1_km_s": dv1,
                "delta_v2_km_s": dv2,
                "delta_v3_km_s": dv3,
                "total_delta_v_km_s": total,
                "hohmann_total_delta_v_km_s": hohmann_total,
                "recommended": "bielliptic" if total < hohmann_total else "hohmann",
                "savings_km_s": hohmann_total - total,
            }
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def combined_plane_change(r_km: float, delta_inc_deg: float,
                              sma_target_km: float = None) -> dict:
        """
        Delta-v of a plane change for a circular orbit of radius r_km. If
        sma_target_km is provided, also computes the cost of a COMBINED
        maneuver (plane change + velocity change at the same point) and
        compares it to the split sequence.

        Args:
            r_km: Orbital radius at the maneuver point (km).
            delta_inc_deg: Desired inclination change (degrees).
            sma_target_km: Semi-major axis of the targeted transfer orbit (km),
                if the maneuver also changes the energy.
        """
        try:
            mu = EARTH_MU_KM3_S2
            di = math.radians(delta_inc_deg)
            v1 = math.sqrt(mu / r_km)
            pure_plane_dv = 2.0 * v1 * math.sin(di / 2.0)
            result = {
                "orbital_velocity_km_s": v1,
                "pure_plane_change_dv_km_s": pure_plane_dv,
            }
            if sma_target_km is not None:
                v2 = math.sqrt(mu * (2.0 / r_km - 1.0 / sma_target_km))
                combined = math.sqrt(v1**2 + v2**2 - 2.0 * v1 * v2 * math.cos(di))
                sequential = abs(v2 - v1) + pure_plane_dv
                result.update({
                    "transfer_velocity_km_s": v2,
                    "combined_dv_km_s": combined,
                    "sequential_dv_km_s": sequential,
                    "recommended": "combined" if combined < sequential else "sequential",
                    "savings_km_s": sequential - combined,
                })
            return result
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def solve_lambert(r1_km: list, r2_km: list, tof_sec: float,
                      mu_km3_s2: float = EARTH_MU_KM3_S2,
                      prograde: bool = True) -> dict:
        """
        Solves Lambert's problem (universal variables): finds the transfer
        orbit connecting two position vectors in a given time of flight.
        Returns the velocity vectors at departure (v1) and arrival (v2).

        Args:
            r1_km: Departure position vector [x, y, z] in km (inertial frame).
            r2_km: Arrival position vector [x, y, z] in km.
            tof_sec: Time of flight in seconds.
            mu_km3_s2: Gravitational parameter of the central body (default Earth).
            prograde: True for a prograde transfer (default), False retrograde.
        """
        try:
            mu = mu_km3_s2
            r1 = list(map(float, r1_km))
            r2 = list(map(float, r2_km))
            r1n = math.sqrt(sum(c*c for c in r1))
            r2n = math.sqrt(sum(c*c for c in r2))
            cx = r1[1]*r2[2] - r1[2]*r2[1]
            cy = r1[2]*r2[0] - r1[0]*r2[2]
            cz = r1[0]*r2[1] - r1[1]*r2[0]
            cos_dnu = max(-1.0, min(1.0, (r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])/(r1n*r2n)))
            dnu = math.acos(cos_dnu)
            if prograde:
                if cz < 0:
                    dnu = 2*math.pi - dnu
            else:
                if cz >= 0:
                    dnu = 2*math.pi - dnu

            A = math.sin(dnu) * math.sqrt(r1n * r2n / (1 - math.cos(dnu)))
            if abs(A) < 1e-12:
                return {"error": "Degenerate geometry (A≈0): Lambert unsolvable."}

            psi = 0.0
            c2, c3 = 0.5, 1.0/6.0
            psi_up, psi_low = 4*math.pi**2, -4*math.pi
            tof_calc = 0.0
            y = r1n + r2n
            for _ in range(1000):
                y = r1n + r2n + A * (psi*c3 - 1) / math.sqrt(c2)
                if A > 0 and y < 0:
                    # raise psi_low until y >= 0
                    while y < 0:
                        psi_low += 0.1
                        psi = (psi_up + psi_low) / 2 if False else psi_low
                        c2, c3 = _stumpff(psi)
                        y = r1n + r2n + A * (psi*c3 - 1) / math.sqrt(c2)
                chi = math.sqrt(y / c2)
                tof_calc = (chi**3 * c3 + A * math.sqrt(y)) / math.sqrt(mu)
                if abs(tof_calc - tof_sec) < 1e-6:
                    break
                if tof_calc <= tof_sec:
                    psi_low = psi
                else:
                    psi_up = psi
                psi = (psi_up + psi_low) / 2.0
                c2, c3 = _stumpff(psi)

            f = 1 - y / r1n
            g = A * math.sqrt(y / mu)
            gdot = 1 - y / r2n
            v1 = [(r2[i] - f*r1[i]) / g for i in range(3)]
            v2 = [(gdot*r2[i] - r1[i]) / g for i in range(3)]
            return {
                "v1_km_s": v1,
                "v2_km_s": v2,
                "v1_mag_km_s": math.sqrt(sum(c*c for c in v1)),
                "v2_mag_km_s": math.sqrt(sum(c*c for c in v2)),
                "transfer_angle_deg": math.degrees(dnu),
                "tof_achieved_sec": tof_calc,
            }
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def time_to_apsis(propagator_name: str, apsis: str = "periapsis") -> dict:
        """
        Time remaining (from the propagator's current state) until the next
        passage at periapsis, apoapsis, ascending node or descending node.
        Useful to position a maneuver at the right place.

        Args:
            propagator_name: Name of the propagator (already propagated at least once).
            apsis: "periapsis" (default), "apoapsis", "ascending_node" or
                "descending_node".
        """
        try:
            return session.time_to_apsis(propagator_name, apsis)
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def convert_time(epoch_utc: str, to_scale: str = "TT") -> dict:
        """
        Converts a UTC epoch (GMAT format "01 Jan 2026 00:00:00.000") to
        another time scale: TAI or TT.

        Offsets: TAI = UTC + ΔAT (37 s since 2017), TT = TAI + 32.184 s.

        Args:
            epoch_utc: UTC epoch in UTCGregorian format.
            to_scale: "TT" (default) or "TAI".
        """
        try:
            from gmat_session import _epoch_to_jd, _jd_to_gregorian
            scale = to_scale.upper()
            if scale not in ("TT", "TAI"):
                return {"error": "to_scale must be 'TT' or 'TAI'."}
            jd_utc = _epoch_to_jd(epoch_utc)
            delta_at = 37.0  # leap seconds accumulated since 2017 (constant since)
            offset = delta_at if scale == "TAI" else delta_at + 32.184
            jd_out = jd_utc + offset / 86400.0
            return {
                "epoch_utc": epoch_utc,
                "scale": scale,
                "offset_sec": offset,
                "epoch_converted": _jd_to_gregorian(jd_out),
                "jd_converted": jd_out,
            }
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def keplerian_to_equinoctial(sma_km: float, ecc: float, inc_deg: float,
                                 raan_deg: float, aop_deg: float,
                                 ta_deg: float) -> dict:
        """
        Converts Keplerian elements to equinoctial elements (a, h, k, p, q,
        lambda), which are more robust near singularities (e≈0 or i≈0).
        """
        try:
            inc = math.radians(inc_deg)
            raan = math.radians(raan_deg)
            aop = math.radians(aop_deg)
            ta = math.radians(ta_deg)
            # mean anomaly
            E = math.atan2(math.sqrt(1-ecc**2)*math.sin(ta), ecc + math.cos(ta))
            M = E - ecc*math.sin(E)
            h = ecc * math.sin(aop + raan)
            k = ecc * math.cos(aop + raan)
            p = math.tan(inc/2.0) * math.sin(raan)
            q = math.tan(inc/2.0) * math.cos(raan)
            lam = math.degrees((raan + aop + M) % (2*math.pi))
            return {"sma_km": sma_km, "h": h, "k": k, "p": p, "q": q,
                    "mean_lon_deg": lam}
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def equinoctial_to_keplerian(sma_km: float, h: float, k: float, p: float,
                                 q: float, mean_lon_deg: float) -> dict:
        """
        Converts equinoctial elements (a, h, k, p, q, lambda) back to
        classical Keplerian elements.
        """
        try:
            ecc = math.sqrt(h**2 + k**2)
            inc = 2.0 * math.atan(math.sqrt(p**2 + q**2))
            raan = math.atan2(p, q)
            aop_raan = math.atan2(h, k)
            aop = aop_raan - raan
            M = math.radians(mean_lon_deg) - aop_raan
            # Kepler: M -> E -> ta
            E = M
            for _ in range(50):
                E = E - (E - ecc*math.sin(E) - M) / (1 - ecc*math.cos(E))
            ta = math.atan2(math.sqrt(1-ecc**2)*math.sin(E), math.cos(E) - ecc)
            return {
                "sma_km": sma_km,
                "ecc": ecc,
                "inc_deg": math.degrees(inc),
                "raan_deg": math.degrees(raan) % 360.0,
                "aop_deg": math.degrees(aop) % 360.0,
                "ta_deg": math.degrees(ta) % 360.0,
            }
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def tsiolkovsky(delta_v_km_s: float = None, isp_s: float = 300.0,
                    mass_initial_kg: float = None, mass_final_kg: float = None,
                    propellant_kg: float = None) -> dict:
        """
        Tsiolkovsky equation. Solves the missing quantity from the others
        (delta-v, masses, propellant) for a given Isp. WARNING: the delta-v
        is in km/s on input but converted to m/s in the exponential.

        Args:
            delta_v_km_s: Delta-v in km/s (provided, or to compute if None).
            isp_s: Specific impulse in seconds (default 300 s).
            mass_initial_kg: Initial (wet) mass in kg.
            mass_final_kg: Final (dry) mass in kg.
            propellant_kg: Propellant mass in kg.
        """
        try:
            ve = isp_s * G0  # m/s
            mi, mf, mp = mass_initial_kg, mass_final_kg, propellant_kg
            if mi is not None and mp is not None and mf is None:
                mf = mi - mp
            elif mf is not None and mp is not None and mi is None:
                mi = mf + mp

            if delta_v_km_s is None:
                if mi is None or mf is None:
                    return {"error": "Provide two of the three masses to compute the delta-v."}
                dv = ve * math.log(mi / mf) / 1000.0
                return {"delta_v_km_s": dv, "delta_v_m_s": dv*1000,
                        "mass_initial_kg": mi, "mass_final_kg": mf,
                        "propellant_kg": mi - mf, "isp_s": isp_s}
            else:
                dv_ms = delta_v_km_s * 1000.0
                ratio = math.exp(dv_ms / ve)  # mi/mf
                if mi is not None:
                    mf = mi / ratio
                elif mf is not None:
                    mi = mf * ratio
                else:
                    return {"error": "Provide at least one mass (initial or final)."}
                return {"delta_v_km_s": delta_v_km_s, "mass_ratio": ratio,
                        "mass_initial_kg": mi, "mass_final_kg": mf,
                        "propellant_kg": mi - mf, "isp_s": isp_s}
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def get_gravitational_parameters() -> dict:
        """
        Table of gravitational parameters μ (km³/s²) and equatorial radii
        (km) of the main bodies of the solar system, to choose the central
        body or size an interplanetary transfer.
        """
        return {
            "mu_km3_s2": dict(MU),
            "radius_km": dict(BODY_RADIUS_KM),
            "g0_m_s2": G0,
        }
