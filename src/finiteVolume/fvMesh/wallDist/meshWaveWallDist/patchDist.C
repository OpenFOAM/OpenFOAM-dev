/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    mesh_(mesh),
    patchIDs_(patchIDs),
    correctWalls_(correctWalls),
    nUnset_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::patchDist::y() const
{
    tmp<volScalarField> ty
    (
        new volScalarField
        (
            IOobject
            (
                "y",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("y", dimLength, GREAT)
        )
    );

    volScalarField& y = ty();

    // Calculate distance starting from patch faces
    patchWave wave(mesh_, patchIDs_, correctWalls_);

    // Transfer cell values from wave into *this
    y.transfer(wave.distance());

    // Transfer values on patches into boundaryField of *this
    forAll(y.boundaryField(), patchI)
    {
        if (!isA<emptyFvPatchScalarField>(y.boundaryField()[patchI]))
        {
            scalarField& waveFld = wave.patchDistance()[patchI];

            y.boundaryField()[patchI].transfer(waveFld);
        }
    }

    // Transfer number of unset values
    nUnset_ = wave.nUnset();

    return ty;
}


// ************************************************************************* //
