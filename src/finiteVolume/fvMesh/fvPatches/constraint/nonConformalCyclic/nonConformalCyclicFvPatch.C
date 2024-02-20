/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "nonConformalCyclicFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, nonConformalCyclicFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicFvPatch::nonConformalCyclicFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    cyclicFvPatch(patch, bm),
    nonConformalCoupledFvPatch(static_cast<const fvPatch&>(*this)),
    nonConformalCyclicPolyPatch_
    (
        refCast<const nonConformalCyclicPolyPatch>(patch)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicFvPatch::~nonConformalCyclicFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalCyclicPolyPatch&
Foam::nonConformalCyclicFvPatch::nonConformalCyclicPatch() const
{
    return nonConformalCyclicPolyPatch_;
}


const Foam::nonConformalCyclicFvPatch&
Foam::nonConformalCyclicFvPatch::nbrPatch() const
{
    return refCast<const nonConformalCyclicFvPatch>(cyclicFvPatch::nbrPatch());
}


Foam::label Foam::nonConformalCyclicFvPatch::start() const
{
    return nonConformalFvPatch::start();
}


Foam::label Foam::nonConformalCyclicFvPatch::size() const
{
    return nonConformalFvPatch::size();
}


bool Foam::nonConformalCyclicFvPatch::coupled() const
{
    return true;
}


const Foam::labelUList& Foam::nonConformalCyclicFvPatch::faceCells() const
{
    return nonConformalFvPatch::faceCells();
}


void Foam::nonConformalCyclicFvPatch::makeWeights(scalarField& w) const
{
    nonConformalCoupledFvPatch::makeWeights
    (
        w,
        nbrPatch().Sf(),
        nbrPatch().coupledFvPatch::delta()
    );
}


Foam::tmp<Foam::vectorField>
Foam::nonConformalCyclicFvPatch::delta() const
{
    return coupledFvPatch::delta(nbrPatch().coupledFvPatch::delta());
}


// ************************************************************************* //
