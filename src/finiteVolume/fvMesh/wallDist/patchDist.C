/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "patchDist.H"
#include "patchWave.H"
#include "fvMesh.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDist::patchDist
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    volScalarField
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    patchIDs_(patchIDs),
    correctWalls_(correctWalls),
    nUnset_(0)
{
    patchDist::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchDist::~patchDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchDist::correct()
{
    // Calculate distance starting from patch faces
    patchWave wave(mesh(), patchIDs_, correctWalls_);

    // Transfer cell values from wave into *this
    transfer(wave.distance());

    // Transfer values on patches into boundaryField of *this
    forAll(boundaryField(), patchI)
    {
        if (!isA<emptyFvPatchScalarField>(boundaryField()[patchI]))
        {
            scalarField& waveFld = wave.patchDistance()[patchI];

            boundaryField()[patchI].transfer(waveFld);
        }
    }

    // Transfer number of unset values
    nUnset_ = wave.nUnset();
}


// ************************************************************************* //
