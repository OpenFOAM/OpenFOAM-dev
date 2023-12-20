/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "fvPatchDistWave.H"
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
    nCorrectors_(dict.lookupOrDefault<label>("nCorrectors", 2)),
    minFaceFraction_(dict.lookupOrDefault<scalar>("minFaceFraction", 1e-1))
{}


Foam::patchDistMethods::meshWave::meshWave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const label nCorrectors,
    const scalar minFaceFraction
)
:
    patchDistMethod(mesh, patchIDs),
    nCorrectors_(nCorrectors),
    minFaceFraction_(minFaceFraction)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::meshWave::correct(volScalarField& y)
{
    y = dimensionedScalar(dimLength, great);

    const label nUnset =
        fvPatchDistWave::calculateAndCorrect
        (
            mesh_,
            patchIndices_,
            minFaceFraction_,
            nCorrectors_,
            y
        );

    // Update coupled and transform BCs
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
        fvPatchDistWave::calculateAndCorrect
        (
            mesh_,
            patchIndices_,
            minFaceFraction_,
            nCorrectors_,
            y,
            n
        );

    // Update coupled and transform BCs
    y.correctBoundaryConditions();
    n.correctBoundaryConditions();

    return nUnset > 0;
}


// ************************************************************************* //
