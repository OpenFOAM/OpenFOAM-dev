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

#include "remove.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(remove, 0);
        addToRunTimeSelectionTable(zoneGenerator, remove, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::remove::remove
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    pointZoneNames_
    (
        dict.found("pointZone")
      ? wordList(1, dict.lookup<word>("pointZone"))
      : dict.lookupOrDefault("pointZones", wordList())
    ),
    cellZoneNames_
    (
        dict.found("cellZone")
      ? wordList(1, dict.lookup<word>("cellZone"))
      : dict.lookupOrDefault("cellZones", wordList())
    ),
    faceZoneNames_
    (
        dict.found("faceZone")
      ? wordList(1, dict.lookup<word>("faceZone"))
      : dict.lookupOrDefault("faceZones", wordList())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::remove::~remove()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::remove::generate() const
{
    forAll(pointZoneNames_, i)
    {
        if (mesh_.pointZones().found(pointZoneNames_[i]))
        {
            const_cast<pointZoneList&>(mesh_.pointZones()).remove
            (
                pointZoneNames_[i]
            );
        }
    }

    forAll(cellZoneNames_, i)
    {
        if (mesh_.cellZones().found(cellZoneNames_[i]))
        {
            const_cast<cellZoneList&>(mesh_.cellZones()).remove
            (
                cellZoneNames_[i]
            );
        }
    }

    forAll(faceZoneNames_, i)
    {
        if (mesh_.faceZones().found(faceZoneNames_[i]))
        {
            const_cast<faceZoneList&>(mesh_.faceZones()).remove
            (
                faceZoneNames_[i]
            );
        }
    }

    return zoneSet();
}


// ************************************************************************* //
