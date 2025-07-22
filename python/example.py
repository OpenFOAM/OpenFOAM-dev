from moldfoam import config

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