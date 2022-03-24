/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "patchDistWave.H"
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
    correctWalls_(dict.lookupOrDefault<Switch>("correctWalls", true))
{}


Foam::patchDistMethods::meshWave::meshWave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(correctWalls)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::meshWave::correct(volScalarField& y)
{
    y = dimensionedScalar(dimLength, great);

    const label nUnset =
        patchDistWave::wave<wallPoint>
        (
            mesh_,
            patchIDs_,
            y.primitiveFieldRef(),
            correctWalls_
        );

    y.correctBoundaryConditions();

    return nUnset > 0;
}


bool Foam::patchDistMethods::meshWave::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    y = dimensionedScalar(dimLength, great);

    const label nUnset =
        patchDistWave::wave<wallPointData<vector>, fvPatchField>
        (
            mesh_,
            patchIDs_,
            n.boundaryField(),
            y.primitiveFieldRef(),
            y.boundaryFieldRef(),
            n.primitiveFieldRef(),
            n.boundaryFieldRef(),
            correctWalls_
        );

    // Update coupled and transform BCs
    y.correctBoundaryConditions();
    n.correctBoundaryConditions();

    return nUnset > 0;
}


// ************************************************************************* //
