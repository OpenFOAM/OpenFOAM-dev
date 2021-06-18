#include "label.H"
#include "HashSet.H"
#include "labelVector.H"
#include "setWriter.H"

namespace Foam
{
    class polyMesh;
    class wedgePolyPatch;
    class surfaceWriter;

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

    //- Write out the weights-sums on all the AMI patches
    void writeAMIWeightsSums(const polyMesh&);

    //- Write out the weights-sum on the given AMI patch
    void writeAMIWeightsSum
    (
        const polyMesh&,
        const primitivePatch&,
        const scalarField&,
        const fileName&
    );

    label checkGeometry
    (
        const polyMesh& mesh,
        const bool allGeometry,
        const autoPtr<surfaceWriter>&,
        const autoPtr<setWriter<scalar>>&
    );
}
