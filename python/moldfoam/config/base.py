from typing import Any, Literal
from pydantic import BaseModel, Field

from moldfoam import j2_env

class FOAMConfig(BaseModel):
    """Base class for OpenFOAM configurations."""

    class Config:
        populate_by_name = True

    foam_version: str = Field(default="2.0", description="FoamFile format version")
    format: Literal['ascii', 'binary'] = Field(default="ascii", description="File format")

    # These vary by file type and should be set by subclasses
    class_: str = Field(..., alias="class", description="OpenFOAM class type")
    location: str = Field(..., description="File location (e.g., 'system', 'constant')")
    object: str = Field(..., description="Object name (usually filename)")

    @property
    def template_name(self) -> str:
        """Return the template name for rendering."""
        raise NotImplementedError("Subclasses must define a template_name")

    def to_openfoam_dict(self) -> dict[str, Any]:
        """Get dictionary with camelCase keys for OpenFOAM templates."""
        return self.dict(by_alias=True)

    def render(self) -> str:
        """Render the controlDict template with the current configuration."""
        template = j2_env.get_template(self.template_name)
        return template.render(**self.to_openfoam_dict())
    
