#include "label.H"
#include "autoPtr.H"

namespace Foam
{
    class polyMesh;
    class surfaceWriter;

    label checkTopology
    (
        const polyMesh&,
        const bool,
        const bool,
        const autoPtr<surfaceWriter>&
    );
}
