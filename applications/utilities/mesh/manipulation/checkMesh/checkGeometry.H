#include "label.H"
#include "HashSet.H"
#include "labelVector.H"

namespace Foam
{
    class polyMesh;
    class wedgePolyPatch;

    label findOppositeWedge(const polyMesh&, const wedgePolyPatch&);

    //- Check wedge orientation
    bool checkWedges
    (
        const polyMesh&,
        const bool report,
        const Vector<label>&,
        labelHashSet*
    );

    //- Check 0th vertex on coupled faces
    bool checkCoupledPoints(const polyMesh&, const bool report, labelHashSet*);

    label checkGeometry(const polyMesh& mesh, const bool allGeometry);
}
