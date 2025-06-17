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

#include "write.H"
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(write, 0);
        addToRunTimeSelectionTable(zoneGenerator, write, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneListType>
void Foam::zoneGenerators::write::writeZoneType(const ZoneListType& zones) const
{
    if
    (
        mesh_.time().timeIndex() == mesh_.time().startTimeIndex()
     || (mesh_.time().writeTime() && zones.timeIndex() > timeIndex_)
    )
    {
        if (mesh_.time().timeIndex() != mesh_.time().startTimeIndex())
        {
            zones.updateTimeInstance();
        }

        zones.write();

        timeIndex_ = mesh_.time().timeIndex();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::write::write
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneType_
    (
        zoneTypesAllNames.lookupOrDefault("zoneType", dict, zoneTypesAll::all)
    ),
    timeIndex_(mesh.time().timeIndex())
{
    moveUpdate_ = dict.lookupOrDefault("moveUpdate", true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::write::~write()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::write::generate() const
{
    if (zoneType_ == zoneTypesAll::point || zoneType_ == zoneTypesAll::all)
    {
        writeZoneType(mesh_.pointZones());
    }

    if (zoneType_ == zoneTypesAll::cell || zoneType_ == zoneTypesAll::all)
    {
        writeZoneType(mesh_.cellZones());
    }

    if (zoneType_ == zoneTypesAll::face || zoneType_ == zoneTypesAll::all)
    {
        writeZoneType(mesh_.faceZones());
    }

    return zoneSet();
}


// ************************************************************************* //
