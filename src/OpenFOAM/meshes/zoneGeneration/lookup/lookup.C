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
            const pointZone* pZonePtr =
                mesh_.pointZones().lookupPtr(zoneName_);

            if (pZonePtr == nullptr)
            {
                FatalIOErrorInFunction(dict_)
                    << "Cannot find pointZone " << zoneName_ << nl
                    << "    Available pointZones "
                    << mesh_.pointZones().sortedToc()
                    << exit(FatalIOError);
            }

            const pointZone& pZone = *pZonePtr;
            moveUpdate_ = moveUpdate_ || pZone.moveUpdate();
            return zoneSet(pZone);
        }

        case zoneTypesAll::cell:
        {
            const cellZone* cZonePtr =
                mesh_.cellZones().lookupPtr(zoneName_);

            if (cZonePtr == nullptr)
            {
                FatalIOErrorInFunction(dict_)
                    << "Cannot find cellZone " << zoneName_ << nl
                    << "    Available cellZones "
                    << mesh_.cellZones().sortedToc()
                    << exit(FatalIOError);
            }

            const cellZone& cZone = *cZonePtr;
            moveUpdate_ = moveUpdate_ || cZone.moveUpdate();
            return zoneSet(cZone);
        }

        case zoneTypesAll::face:
        {
            const faceZone* fZonePtr =
                mesh_.faceZones().lookupPtr(zoneName_);

            if (fZonePtr == nullptr)
            {
                FatalIOErrorInFunction(dict_)
                    << "Cannot find faceZone " << zoneName_ << nl
                    << "    Available faceZones "
                    << mesh_.faceZones().sortedToc()
                    << exit(FatalIOError);
            }

            const faceZone& fZone = *fZonePtr;
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
                zs.pZone = mesh_.pointZones()[zoneName_];
                found = true;
                moveUpdate_ = moveUpdate_ || zs.pZone().moveUpdate();
            }

            if (mesh_.cellZones().found(zoneName_))
            {
                zs.cZone = mesh_.cellZones()[zoneName_];
                found = true;
                moveUpdate_ = moveUpdate_ || zs.cZone().moveUpdate();
            }

            if (mesh_.faceZones().found(zoneName_))
            {
                zs.fZone = mesh_.faceZones()[zoneName_];
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
