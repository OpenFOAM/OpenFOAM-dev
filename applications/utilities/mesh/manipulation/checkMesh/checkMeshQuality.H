#include "polyMesh.H"

namespace Foam
{
    class surfaceWriter;

    label checkMeshQuality
    (
        const polyMesh&,
        const dictionary&,
        const autoPtr<surfaceWriter>&
    );
}
