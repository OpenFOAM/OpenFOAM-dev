/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nonConformalCoupledFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "nonConformalPolyFacesFvsPatchLabelField.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCoupledFvPatch, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nonConformalCoupledFvPatch::makeWeights
(
    scalarField& w,
    const vectorField& nbrSf,
    const vectorField& nbrDelta
) const
{
    const vectorField delta(patch_.coupledFvPatch::delta());

    const scalarField nfDelta(patch_.nf() & delta);

    const scalarField nbrNfDelta((nbrSf/(mag(nbrSf) + vSmall)) & nbrDelta);

    forAll(delta, facei)
    {
        const scalar ndoi = nfDelta[facei];
        const scalar ndni = nbrNfDelta[facei];
        const scalar ndi = ndoi + ndni;

        // !!! Note this is a different form of stabilisation compared to
        // coupledFvPatch. This is necessary to prevent negative weights on
        // very small faces. Should this also be the standard form of
        // stabilisation in coupledFvPatch and in surfaceInterpolation?

        if (ndoi > vSmall && ndni > vSmall)
        {
            w[facei] = ndni/ndi;
        }
        else
        {
            const scalar doi = mag(delta[facei]);
            const scalar dni = mag(nbrDelta[facei]);
            const scalar di = doi + dni;

            w[facei] = dni/di;
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledFvPatch::nonConformalCoupledFvPatch
(
    const fvPatch& patch
)
:
    nonConformalFvPatch(patch),
    patch_(refCast<const coupledFvPatch>(patch)),
    nonConformalCoupledPolyPatch_
    (
        refCast<const nonConformalCoupledPolyPatch>(patch.patch())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledFvPatch::~nonConformalCoupledFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalCoupledPolyPatch&
Foam::nonConformalCoupledFvPatch::nonConformalCoupledPatch() const
{
    return nonConformalCoupledPolyPatch_;
}


bool Foam::nonConformalCoupledFvPatch::owner() const
{
    return nonConformalCoupledPolyPatch_.owner();
}


bool Foam::nonConformalCoupledFvPatch::neighbour() const
{
    return nonConformalCoupledPolyPatch_.neighbour();
}


const Foam::transformer& Foam::nonConformalCoupledFvPatch::transform() const
{
    return nonConformalCoupledPolyPatch_.transform();
}


const Foam::word& Foam::nonConformalCoupledFvPatch::errorPatchName() const
{
    return nonConformalCoupledPolyPatch_.errorPatchName();
}


Foam::label Foam::nonConformalCoupledFvPatch::errorPatchIndex() const
{
    return nonConformalCoupledPolyPatch_.errorPatchIndex();
}


const Foam::nonConformalErrorFvPatch&
Foam::nonConformalCoupledFvPatch::errorPatch() const
{
    return
        refCast<const nonConformalErrorFvPatch>
        (
            patch_.boundaryMesh()[errorPatchIndex()]
        );
}


Foam::word Foam::nonConformalCoupledFvPatch::polyFacesType() const
{
    return nonConformalPolyFacesFvsPatchLabelField::typeName;
}


// ************************************************************************* //
