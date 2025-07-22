from pydantic import Field, validator, model_validator
from typing import Literal, Optional

from ..base import FOAMConfig

class ControlConfig(FOAMConfig):
    """OpenFOAM controlDict configuration."""

    @property
    def template_name(self) -> str:
        return "system/controlDict.jinja"
    
    # FoamFile header
    class_: str = Field(default="dictionary")
    location: str = Field(default="system")
    object: str = Field(default="controlDict")

    # Required parameters
    solver: str = Field(..., description="OpenFOAM solver to use")
    end_time: float = Field(..., alias="endTime", description="Simulation end time")
    delta_t: float = Field(..., alias="deltaT", description="Time step size")
    write_interval: float = Field(..., alias="writeInterval", description="Write interval")
    adjust_time_step: bool = Field(..., alias="adjustTimeStep", description="Adjust time step")

    # Conditionally required fields (when adjust_time_step is True)
    max_co: Optional[float] = Field(default=1.0, alias="maxCo", description="Maximum Courant number")
    max_alpha_co: Optional[float] = Field(default=1.0, alias="maxAlphaCo", description="Maximum interface Courant")
    max_delta_t: Optional[float] = Field(default=1.0, alias="maxDeltaT", description="Maximum time step")

    # Optional parameters
    application: str = Field(default='foamRun')
    start_from: Literal['startTime', 'firstTime', 'latestTime'] = Field(
        default='startTime', alias="startFrom"
    )
    start_time: float = Field(default=0, alias="startTime")
    stop_at: Literal['endTime', 'writeTime', 'noWriteNow', 'nextWrite'] = Field(
        default='endTime', alias="stopAt"
    )
    write_control: Literal['timeStep', 'runTime', 'adjustableRunTime', 'cpuTime', 'clockTime'] = Field(
        default='adjustableRunTime', alias="writeControl"
    )
    purge_write: int = Field(default=0, alias="purgeWrite")
    write_format: Literal['ascii', 'binary'] = Field(default='ascii', alias="writeFormat")
    write_precision: int = Field(default=6, alias="writePrecision")
    write_compression: Literal['on', 'off'] = Field(default='off', alias="writeCompression")
    time_format: Literal['fixed', 'scientific', 'general'] = Field(default='general', alias="timeFormat")
    time_precision: int = Field(default=6, alias="timePrecision")
    run_time_modifiable: bool = Field(default=True, alias="runTimeModifiable")
    track_pressure: bool = Field(default=True, alias="trackPressure")

    @validator('delta_t', 'max_delta_t')
    def validate_positive_time(cls, v):
        if v <= 0:
            raise ValueError('Time values must be positive')
        return v

    @model_validator(mode='before')
    def validate_conditional_requirements(cls, values):
        """Validate that required fields are present when adjust_time_step=True."""
        adjust_time_step = values.get('adjust_time_step')
        
        if adjust_time_step:
            required_fields = ['max_co', 'max_alpha_co', 'max_delta_t']
            missing_fields = []
            
            for field in required_fields:
                if field not in values or values[field] is None:
                    missing_fields.append(field)
            
            if missing_fields:
                raise ValueError(
                    f"When adjust_time_step=True, these fields are required: {missing_fields}"
                )
                
        return values