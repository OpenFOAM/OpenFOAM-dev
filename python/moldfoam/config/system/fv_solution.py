from pydantic import BaseModel, Field, validator, model_validator
from typing import Literal, Optional

from ..base import FOAMConfig, FoamConfigBase


class PIMPLEConfig(FoamConfigBase):
    momentum_predictor: Literal['yes', 'no'] = Field(alias="momentumPredictor")
    n_outer_correctors: int = Field(alias="nOuterCorrectors")
    n_correctors: int = Field(alias="nCorrectors")
    n_non_orthogonal_correctors: int = Field(alias="nNonOrthogonalCorrectors")


class PressureConfig(FoamConfigBase):
    preconditioner: Literal["DIC", "GAMG"] = Field()
    atol: float = Field(default=1e-8, alias="tolerance")
    rtol: float = Field(default=1e-2, alias="relTol")
    max_iter: int = Field(default=100, alias="maxIter")


class FVSolutionConfig(FOAMConfig):
    """OpenFOAM fvSolution configuration."""

    @property
    def template_name(self) -> str:
        return "system/fvSolution.jinja"

    # FoamFile header
    class_: str = Field(default="dictionary")
    location: str = Field(default="system")
    object: str = Field(default="fvSolution")

    type_: Literal['VoF', 'default'] = Field()
    pimple: PIMPLEConfig = Field(alias="PIMPLE")
    pressure: PressureConfig = Field(alias="pressure")