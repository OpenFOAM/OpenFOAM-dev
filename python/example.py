from moldfoam import config, dimensions
from moldfoam import field

template = config.system.control.ControlConfig(
    solver="compressibleVoF",
    end_time=0.4,
    delta_t=0.01,
    write_interval=0.01,
    adjust_time_step=False,
    max_co=1000.0,
    max_alpha_co=10.0,
    max_delta_t=0.01,
)

print(template.render())

template = config.system.decompose_par.DecomposeParConfig(
    num_subdomains=8,
)

print(template.render())


template = config.system.fv_schemes.FVSchemesConfig(
    type_='VoF',
    interfaceCompression=0,
)

print(template.render())


template = config.system.fv_solution.FVSolutionConfig(
    type_='VoF',
    pimple=config.system.fv_solution.PIMPLEConfig(
        momentum_predictor='no',
        n_outer_correctors=1,
        n_correctors=2,
        n_non_orthogonal_correctors=1,
    ),
    pressure=config.system.fv_solution.PressureConfig(
        preconditioner='GAMG',
    )
)

print(template.to_openfoam_dict())

print(template.render())

T_wall = 323.15  # Wall temperature [K]
T0 = 500.0  # Initial temperature [K]
h = 1200.0  # Heat flux [W/m^2/K]

# Boundary conditions
bcs = {
    "\"MOLD_PART.*\"": field.UniformBC(
        type="externalWallHeatFluxTemperature",
        value=T_wall,
        additional_vars={
            "Ta": T_wall,  # Ambient temperature
            "h": f"uniform {h}",
        }
    ),
    "MELT_POINT": field.UniformBC(
        type="fixedValue",
        value=T0,
    ),
    "MELT_VENT": field.UniformBC(
        type="inletOutlet",
        value=T0,
        additional_vars={
            "inletValue": T_wall,
        }
    )
}

template = field.FOAMField(
    class_='volScalarField',
    object='T',
    location='0',
    dimensions=dimensions.TEMPERATURE,
    is_uniform=True,
    value=500,
    boundary_fields=bcs,
)

print(template.render())

# Pressure field
p_atm = 1e5
bcs = {
    "\"MOLD_PART.*\"": field.UniformBC(
        type="zeroGradient",
        value=p_atm,
    ),
    "MELT_POINT": field.UniformBC(
        type="zeroGradient",
        value=p_atm,
    ),
    # TODO: Need pressure vent BCs
    # "MELT_VENT": field.UniformBC(
    #     type="codedMixed",
    #     value=T0,
    #     additional_vars={
    #         "inletValue": T_wall,
    #     }
    # )
}

template = field.FOAMField(
    class_='volScalarField',
    object='p_rgh',
    location='0',
    dimensions=dimensions.PRESSURE,
    is_uniform=True,
    value=p_atm,
    boundary_fields=bcs,
)

print(template.render())