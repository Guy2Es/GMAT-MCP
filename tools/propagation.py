# tools/propagation.py
"""
MCP tools related to orbital propagation.
"""

from gmat_session import session


def register_propagation_tools(mcp):

    @mcp.tool()
    def create_propagator(
        spacecraft_name: str,
        step_size_sec: float = 60.0,
        gravity_degree: int = 0,
        gravity_order: int = None,
        gravity_file: str = "JGM3",
        third_bodies: list = None,
        enable_drag: bool = False,
        atmosphere_model: str = "MSISE90",
        enable_srp: bool = False,
        dry_mass_kg: float = 850.0,
        cd: float = 2.2,
        drag_area_m2: float = 15.0,
        cr: float = 1.8,
        srp_area_m2: float = 1.0,
    ) -> dict:
        """
        Creates a propagator (force model + PrinceDormand78 integrator)
        for an already-created spacecraft.

        By default the model is a simple Earth point-mass (2-body),
        consistent with the analytic calculations. Enrich the model according
        to the mission: J2+ harmonics mandatory for SSO/frozen orbits,
        third-bodies for GEO and interplanetary, drag for orbital lifetime in
        low orbit. A common realistic value is gravity_degree=8.

        Args:
            spacecraft_name: Name of the existing spacecraft to propagate.
            step_size_sec: Integration step in seconds (default 60s).
            gravity_degree: Degree of the Earth spherical harmonics.
                0 (default) = simple point-mass; >=2 enables the field (J2 at
                degree=2). Max 70 (JGM2/JGM3) or 360 (EGM96).
            gravity_order: Order of the harmonics (default = gravity_degree).
                Set 0 to keep only the zonal terms (J2, J3...).
            gravity_file: Potential model — "JGM3" (default), "JGM2" or
                "EGM96".
            third_bodies: List of perturbers, e.g. ["Luna", "Sun"]
                ("Moon" accepted as an alias of "Luna").
            enable_drag: Enables atmospheric drag (lifetime / re-entry).
            atmosphere_model: "MSISE90" (default), "JacchiaRoberts" or
                "Exponential". Constant solar/geomagnetic flux (F10.7=150).
            enable_srp: Enables solar radiation pressure (dominant in GEO).
            dry_mass_kg: Spacecraft dry mass (used if drag/SRP).
            cd: Drag coefficient (~2.2 typical).
            drag_area_m2: Drag area in m².
            cr: SRP reflectivity coefficient (~1.8).
            srp_area_m2: Area exposed to SRP in m².
        """
        try:
            return session.create_propagator(
                spacecraft_name,
                step_size_sec=step_size_sec,
                gravity_degree=gravity_degree,
                gravity_order=gravity_order,
                gravity_file=gravity_file,
                third_bodies=third_bodies,
                enable_drag=enable_drag,
                atmosphere_model=atmosphere_model,
                enable_srp=enable_srp,
                dry_mass_kg=dry_mass_kg,
                cd=cd,
                drag_area_m2=drag_area_m2,
                cr=cr,
                srp_area_m2=srp_area_m2,
            )
        except Exception as e:
            return {"error": str(e)}

    @mcp.tool()
    def propagate(propagator_name: str, duration_sec: float) -> dict:
        """
        Propagates an orbit for a given duration and returns the
        spacecraft's final orbital state.

        Args:
            propagator_name: Name of the propagator (created via create_propagator).
            duration_sec: Propagation duration in seconds.
        """
        try:
            return session.propagate(propagator_name, duration_sec)
        except Exception as e:
            return {"error": str(e)}