/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "nearPatchCells_zoneGenerator.H"
#include "patchDistWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(nearPatchCells, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            nearPatchCells,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::nearPatchCells::nearPatchCells
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    patchSet_(mesh.boundaryMesh().patchSet(dict)),
    distance_(dict.lookup<scalar>("distance", dimLength))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::nearPatchCells::~nearPatchCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::nearPatchCells::generate() const
{
    scalarField cellWallDistance(mesh_.nCells());

    patchDistWave::calculate(mesh_, patchSet_, cellWallDistance);

    return zoneSet
    (
        new cellZone
        (
            zoneName_,
            selectIndices
            (
                cellWallDistance,
                [&](const scalar d)
                {
                    return d <= distance_;
                }
            ),
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
