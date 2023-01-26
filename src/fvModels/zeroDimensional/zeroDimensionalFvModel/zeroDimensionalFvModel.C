/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "zeroDimensionalFvModel.H"
#include "fluidThermo.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(zeroDimensionalFvModel, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalFvModel::zeroDimensionalFvModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict)
{
    if (mesh.nGeometricD() != 0)
    {
        FatalIOErrorInFunction(dict)
            << "Zero-dimensional fvModel applied to a "
            << mesh.nGeometricD() << "-dimensional mesh"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalFvModel::~zeroDimensionalFvModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::zeroDimensionalFvModel::movePoints()
{
    return true;
}


void Foam::fv::zeroDimensionalFvModel::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fv::zeroDimensionalFvModel::mapMesh
(
    const polyMeshMap& map
)
{}


void Foam::fv::zeroDimensionalFvModel::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
