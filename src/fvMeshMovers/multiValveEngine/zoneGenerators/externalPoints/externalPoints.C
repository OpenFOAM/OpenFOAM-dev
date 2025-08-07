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

#include "externalPoints.H"
#include "multiValveEngine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(externalPoints, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            externalPoints,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::externalPoints::externalPoints
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::externalPoints::~externalPoints()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::externalPoints::generate() const
{
    const fvMeshMovers::multiValveEngine& mve
    (
        refCast<const fvMeshMovers::multiValveEngine>
        (
            refCast<const fvMesh>(mesh_).mover()
        )
    );

    const fvMeshMovers::multiValveEngine::pistonObject& piston = mve.piston;

    return zoneGenerator::New
    (
        name_,
        mesh_,
        dictionary::entries
        (
            "type", "invert",
            "internalPoints", dictionary::entries
            (
                "type", "point",
                "internalCells", dictionary::entries
                (
                    "type", "cylinder",
                    "zoneType", "cell",
                    "point1", piston.centre() - 10*piston.bore()*piston.axis,
                    "point2", piston.centre() + 10*piston.bore()*piston.axis,
                    "radius", piston.bore()/2
                )
            )
        )
    )->generate();
}


// ************************************************************************* //
