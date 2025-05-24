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

#include "set.H"
#include "polyMesh.H"
#include "pointSet.H"
#include "cellSet.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(set, 0);
        addToRunTimeSelectionTable(zoneGenerator, set, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::set::set
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneType_(zoneTypesNames.read(dict.lookup("zoneType"))),
    setName_(dict.lookupOrDefault("setName", zoneName_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::set::~set()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::set::generate() const
{
    switch (zoneType_)
    {
        case pointZoneType:
        {
            return zoneSet
            (
                new pointZone
                (
                    zoneName_,
                    pointSet(mesh_, setName_).toc(),
                    mesh_.pointZones(),
                    false,
                    false
                )
            );
        }

        case cellZoneType:
        {
            return zoneSet
            (
                new cellZone
                (
                    zoneName_,
                    cellSet(mesh_, setName_).toc(),
                    mesh_.cellZones(),
                    false,
                    false
                )
            );
        }

        case faceZoneType:
        {
            labelList faceIndices(faceSet(mesh_, setName_).toc());

            return zoneSet
            (
                new faceZone
                (
                    zoneName_,
                    faceIndices,
                    boolList(faceIndices.size(), false),
                    mesh_.faceZones(),
                    false,
                    false
                )
            );
        }
    }

    return zoneSet();
}


// ************************************************************************* //
