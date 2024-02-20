/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "nonConformalProcessorCyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalProcessorCyclicFvPatch, 0);
    addToRunTimeSelectionTable
    (
        fvPatch,
        nonConformalProcessorCyclicFvPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicFvPatch::nonConformalProcessorCyclicFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    processorCyclicFvPatch(patch, bm),
    nonConformalCoupledFvPatch(static_cast<const fvPatch&>(*this)),
    nonConformalProcessorCyclicPolyPatch_
    (
        refCast<const nonConformalProcessorCyclicPolyPatch>(patch)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicFvPatch::~nonConformalProcessorCyclicFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalProcessorCyclicPolyPatch&
Foam::nonConformalProcessorCyclicFvPatch::
nonConformalProcessorCyclicPatch() const
{
    return nonConformalProcessorCyclicPolyPatch_;
}


Foam::label Foam::nonConformalProcessorCyclicFvPatch::start() const
{
    return nonConformalFvPatch::start();
}


Foam::label Foam::nonConformalProcessorCyclicFvPatch::size() const
{
    return nonConformalFvPatch::size();
}


bool Foam::nonConformalProcessorCyclicFvPatch::coupled() const
{
    return Pstream::parRun();
}


const Foam::labelUList&
Foam::nonConformalProcessorCyclicFvPatch::faceCells() const
{
    return nonConformalFvPatch::faceCells();
}


void Foam::nonConformalProcessorCyclicFvPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        nonConformalCoupledFvPatch::makeWeights
        (
            w,
          - boundaryMesh().mesh().Sf().boundaryField()[index()],
            boundaryMesh().mesh().Cf().boundaryField()[index()]
          - boundaryMesh().mesh().C().boundaryField()[index()]
        );
    }
    else
    {
        w = 1;
    }
}


Foam::tmp<Foam::vectorField>
Foam::nonConformalProcessorCyclicFvPatch::delta() const
{
    if (Pstream::parRun())
    {
        return
            coupledFvPatch::delta
            (
                boundaryMesh().mesh().Cf().boundaryField()[index()]
              - boundaryMesh().mesh().C().boundaryField()[index()]
            );
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


// ************************************************************************* //
