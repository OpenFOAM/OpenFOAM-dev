from pydantic import Field, validator, model_validator
from typing import Literal, Optional

from ..base import FOAMConfig

class FVSchemesConfig(FOAMConfig):
    """OpenFOAM fvSchemes configuration."""

    @property
    def template_name(self) -> str:
        return "system/fvSchemes.jinja"

    # FoamFile header
    class_: str = Field(default="dictionary")
    location: str = Field(default="system")
    object: str = Field(default="fvSchemes")

    # Required parameters
    type_: Literal['VoF', 'default'] = Field(..., description="Type of fvSchemes configuration")

    interface_compression: Optional[float] = Field(
        default=0.0,
        alias="interfaceCompression",
        description="Interface compression factor for VoF schemes",
    )