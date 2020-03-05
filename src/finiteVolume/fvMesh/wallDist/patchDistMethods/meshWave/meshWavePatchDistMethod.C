/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "meshWavePatchDistMethod.H"
#include "fvMesh.H"
#include "volFields.H"
#include "patchWave.H"
#include "patchDataWave.H"
#include "wallPointData.H"
#include "emptyFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(meshWave, 0);
    addToRunTimeSelectionTable(patchDistMethod, meshWave, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::meshWave::meshWave
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(dict.lookupOrDefault<Switch>("correctWalls", true)),
    nUnset_(0)
{}


Foam::patchDistMethods::meshWave::meshWave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(correctWalls),
    nUnset_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::meshWave::correct(volScalarField& y)
{
    y = dimensionedScalar(dimLength, great);

    // Calculate distance starting from patch faces
    patchWave wave(mesh_, patchIDs_, correctWalls_);

    // Transfer cell values from wave into y
    y.transfer(wave.distance());

    // Transfer values on patches into boundaryField of y
    volScalarField::Boundary& ybf = y.boundaryFieldRef();

    forAll(ybf, patchi)
    {
        if (!isA<emptyFvPatchScalarField>(ybf[patchi]))
        {
            scalarField& waveFld = wave.patchDistance()[patchi];

            ybf[patchi].transfer(waveFld);
        }
    }

    // Transfer number of unset values
    nUnset_ = wave.nUnset();

    return nUnset_ > 0;
}


bool Foam::patchDistMethods::meshWave::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    y = dimensionedScalar(dimLength, great);

    // Collect pointers to data on patches
    UPtrList<vectorField> patchData(mesh_.boundaryMesh().size());

    volVectorField::Boundary& nbf = n.boundaryFieldRef();

    forAll(nbf, patchi)
    {
        patchData.set(patchi, &nbf[patchi]);
    }

    // Do mesh wave
    patchDataWave<wallPointData<vector>> wave
    (
        mesh_,
        patchIDs_,
        patchData,
        correctWalls_
    );

    // Transfer cell values from wave into y and n
    y.transfer(wave.distance());

    n.transfer(wave.cellData());

    // Transfer values on patches into boundaryField of y and n
    volScalarField::Boundary& ybf = y.boundaryFieldRef();

    forAll(ybf, patchi)
    {
        scalarField& waveFld = wave.patchDistance()[patchi];

        if (!isA<emptyFvPatchScalarField>(ybf[patchi]))
        {
            ybf[patchi].transfer(waveFld);

            vectorField& wavePatchData = wave.patchData()[patchi];

            nbf[patchi].transfer(wavePatchData);
        }
    }

    // Update coupled BCs
    y.correctBoundaryConditions();

    // Transfer number of unset values
    nUnset_ = wave.nUnset();

    return nUnset_ > 0;
}


// ************************************************************************* //
