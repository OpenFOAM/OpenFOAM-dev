/*======================================================================*/

// This routine evaluates q(k) (McComb, p61) by summing all wavevectors
// in a k-shell. Then we divide through by the volume of the box
// - to be accurate, by the volume of the ellipsoid enclosing the
// box (assume cells even in each direction). Finally, multiply the
// q(k) values by k^2 to give the full power spectrum E(k). Integrating
// this over the whole range gives the energy in turbulence.

/*======================================================================*/

#include "calcEk.H"
#include "fft.H"
#include "Kmesh.H"
#include "kShellIntegration.H"
#include "volFields.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

graph calcEk
(
    const volVectorField& U,
    const Kmesh& K
)
{
    return kShellIntegration
    (
        fft::forwardTransform
        (
            ReComplexField(U.primitiveField()),
            K.nn()
        ),
        K
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
