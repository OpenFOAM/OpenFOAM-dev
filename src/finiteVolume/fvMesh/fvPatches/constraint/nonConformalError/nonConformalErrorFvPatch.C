/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "nonConformalErrorFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalErrorFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, nonConformalErrorFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalErrorFvPatch::nonConformalErrorFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    fvPatch(patch, bm),
    nonConformalFvPatch(static_cast<const fvPatch&>(*this))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalErrorFvPatch::~nonConformalErrorFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::nonConformalErrorFvPatch::polyFaces() const
{
    return nonConformalFvPatch::polyFaces();
}


Foam::label Foam::nonConformalErrorFvPatch::start() const
{
    return nonConformalFvPatch::start();
}


Foam::label Foam::nonConformalErrorFvPatch::size() const
{
    return nonConformalFvPatch::size();
}


const Foam::labelUList& Foam::nonConformalErrorFvPatch::faceCells() const
{
    return nonConformalFvPatch::faceCells();
}


Foam::tmp<Foam::vectorField> Foam::nonConformalErrorFvPatch::delta() const
{
    return Cf() - Cn();
}


// ************************************************************************* //
