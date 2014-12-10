#include "writePatchGraph.H"
#include "volFields.H"
#include "fvMesh.H"
#include "graph.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writePatchGraph
(
    const volScalarField& vsf,
    const label patchLabel,
    const direction d,
    const word& graphFormat
)
{
    graph
    (
        vsf.name(),
        "position",
        vsf.name(),
        vsf.mesh().boundary()[patchLabel].Cf().component(d),
        vsf.boundaryField()[patchLabel]
    ).write(vsf.time().timePath()/vsf.name(), graphFormat);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
