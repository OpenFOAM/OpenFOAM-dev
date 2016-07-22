#include "label.H"
#include "autoPtr.H"
#include "writer.H"

namespace Foam
{
    class polyMesh;
    class surfaceWriter;

    label checkTopology
    (
        const polyMesh&,
        const bool,
        const bool,
        const autoPtr<surfaceWriter>&,
        const autoPtr<writer<scalar>>&
    );
}
