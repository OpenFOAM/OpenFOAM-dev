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

#include "all.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(all, 0);
        addToRunTimeSelectionTable(zoneGenerator, all, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::all::all
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneType_(zoneTypesNames.read(dict.lookup("zoneType")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::all::~all()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::all::generate() const
{
    switch (zoneType_)
    {
        case zoneTypes::point:
        {
            return zoneSet
            (
                new pointZone
                (
                    zoneName_,
                    identityMap(mesh_.nPoints()),
                    mesh_.pointZones(),
                    false,
                    false
                )
            );
        }

        case zoneTypes::cell:
        {
            return zoneSet
            (
                new cellZone
                (
                    zoneName_,
                    identityMap(mesh_.nCells()),
                    mesh_.cellZones(),
                    false,
                    false
                )
            );
        }

        case zoneTypes::face:
        {
            return zoneSet
            (
                new faceZone
                (
                    zoneName_,
                    identityMap(mesh_.nFaces()),
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
