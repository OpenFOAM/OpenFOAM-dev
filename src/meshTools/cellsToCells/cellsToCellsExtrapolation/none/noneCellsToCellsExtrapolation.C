/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "noneCellsToCellsExtrapolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellsExtrapolations
{
    defineTypeNameAndDebug(none, 0);
    addToRunTimeSelectionTable(cellsToCellsExtrapolation, none, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolations::none::none()
:
    cellsToCellsExtrapolation()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolations::none::~none()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellsToCellsExtrapolations::none::update
(
    const polyMesh& mesh,
    const PackedBoolList& cellCoupleds
)
{
    cellsToCellsExtrapolation::update(cellCoupleds);

    if (extrapolation_)
    {
        FatalErrorInFunction
            << "Mapping of the cells of mesh '" << mesh.name() << "' is "
            << "incomplete and requires extrapolation, but extrapolation is "
            << "deactivated because the selected extrapolation engine is of "
            << "type '" << typeName << "'. Either improve the correspondence "
            << "between the meshes, or select a functional extrapolation "
            << "engine." << exit(FatalError);
    }
}


// ************************************************************************* //
