import moldfoam

from moldfoam.config.control_dict import ControlDictConfig

config = ControlDictConfig(
    solver="compressibleVoF",
    end_time=0.4,
    delta_t=0.01,
    write_interval=0.01,
    adjust_time_step=False,
    max_co=1000.0,
    max_alpha_co=10.0,
    max_delta_t=0.01,
)

print(config)
print(config.to_openfoam_dict())

print(config.render())