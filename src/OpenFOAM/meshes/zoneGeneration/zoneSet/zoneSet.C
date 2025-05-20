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

#include "zoneSet.H"
#include "pointZoneList.H"
#include "cellZoneList.H"
#include "faceZoneList.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneSet::store() const
{
    zoneSet zs;

    if (!pZone.empty())
    {
        pointZone* pZonePtr = pZone.ptr();
        zs.pZone = pZonePtr->zones().append(pZonePtr);
    }
    else if (pZone.valid())
    {
        zs.pZone = pZone();
    }

    if (!cZone.empty())
    {
        cellZone* cZonePtr = cZone.ptr();
        zs.cZone = cZonePtr->zones().append(cZonePtr);
    }
    else if (cZone.valid())
    {
        zs.cZone = cZone();
    }

    if (!fZone.empty())
    {
        faceZone* fZonePtr = fZone.ptr();
        zs.fZone = fZonePtr->zones().append(fZonePtr);
    }
    else if (fZone.valid())
    {
        zs.fZone = fZone();
    }

    return zs;
}


void Foam::zoneSet::operator=(const zoneSet& zs)
{
    if (!zs.pZone.empty())
    {
        pZone = zs.pZone;
    }

    if (!zs.cZone.empty())
    {
        cZone = zs.cZone;
    }

    if (!zs.fZone.empty())
    {
        fZone = zs.fZone;
    }
}


// ************************************************************************* //
