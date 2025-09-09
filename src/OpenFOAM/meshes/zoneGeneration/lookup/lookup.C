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

#include "lookup.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(lookup, 0);
        addToRunTimeSelectionTable(zoneGenerator, lookup, dictionary);
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ZoneListType>
const typename ZoneListType::zoneType& Foam::zoneGenerators::lookup::lookupZone
(
    const ZoneListType& zones,
    const word& zoneName
) const
{
    typedef typename ZoneListType::zoneType ZoneType;

    const ZoneType* zonePtr = zones.lookupPtr(zoneName_);

    if (zonePtr == nullptr)
    {
        FatalIOErrorInFunction(dict_)
            << "Cannot find " << ZoneType::typeName << " " << zoneName_ << nl;

        if (mesh_.pointZones().lookupPtr(zoneName_) != nullptr)
        {
            FatalIOError << "    Found pointZone " << zoneName << nl;
        }

        if (mesh_.cellZones().lookupPtr(zoneName_) != nullptr)
        {
            FatalIOError << "    Found cellZone " << zoneName << nl;
        }

        if (mesh_.faceZones().lookupPtr(zoneName_) != nullptr)
        {
            FatalIOError << "    Found faceZone " << zoneName << nl;
        }

        FatalIOError
            << "    Available " << ZoneType::typeName << "s: "
            << zones.sortedToc()
            << exit(FatalIOError);
    }

    return *zonePtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::lookup::lookup
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::lookup::~lookup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::lookup::generate() const
{
    switch (zoneType_)
    {
        case zoneTypesAll::point:
        {
            const pointZone& pZone = lookupZone(mesh_.pointZones(), zoneName_);
            moveUpdate_ = moveUpdate_ || pZone.moveUpdate();
            return zoneSet(pZone);
        }

        case zoneTypesAll::cell:
        {
            const cellZone& cZone = lookupZone(mesh_.cellZones(), zoneName_);
            moveUpdate_ = moveUpdate_ || cZone.moveUpdate();
            return zoneSet(cZone);
        }

        case zoneTypesAll::face:
        {
            const faceZone& fZone = lookupZone(mesh_.faceZones(), zoneName_);
            moveUpdate_ = moveUpdate_ || fZone.moveUpdate();
            return zoneSet(fZone);
        }

        case zoneTypesAll::all:
        {
            if (zoneName_ == "all")
            {
                FatalIOErrorInFunction(dict_)
                    << "zoneType not specified or set to all for zone all"
                    << exit(FatalIOError);
            }

            zoneSet zs;
            bool found = false;

            if (mesh_.pointZones().found(zoneName_))
            {
                zs = mesh_.pointZones()[zoneName_];
                found = true;
                moveUpdate_ = moveUpdate_ || zs.pZone().moveUpdate();
            }

            if (mesh_.cellZones().found(zoneName_))
            {
                zs = mesh_.cellZones()[zoneName_];
                found = true;
                moveUpdate_ = moveUpdate_ || zs.cZone().moveUpdate();
            }

            if (mesh_.faceZones().found(zoneName_))
            {
                zs = mesh_.faceZones()[zoneName_];
                found = true;
                moveUpdate_ = moveUpdate_ || zs.fZone().moveUpdate();
            }

            if (!found)
            {
                FatalIOErrorInFunction(dict_)
                    << "Cannot find zone " << zoneName_ << exit(FatalIOError);
            }

            return zs;
        }
    }

    return zoneSet();
}


// ************************************************************************* //
