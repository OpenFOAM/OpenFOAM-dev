from pydantic import Field

from ..base import FOAMConfig

class DecomposeParConfig(FOAMConfig):
    """OpenFOAM decomposeParDict configuration."""

    @property
    def template_name(self) -> str:
        return "system/decomposeParDict.jinja"

    # FoamFile header
    class_: str = Field(default="dictionary")
    location: str = Field(default="system")
    object: str = Field(default="decomposeParDict")

    # Required parameters
    num_subdomains: int = Field(..., alias="numberOfSubdomains", description="Number of subdomains for decomposition")