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

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

const Foam::NamedEnum<Foam::zoneTypes, 3>
Foam::zoneTypesNames
{
    "point",
    "cell",
    "face"
};

const Foam::NamedEnum<Foam::zoneTypesAll, 4>
Foam::zoneTypesAllNames
{
    "point",
    "cell",
    "face",
    "all"
};


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneSet::store() const
{
    zoneSet zs;

    if (!pZone_.empty())
    {
        pointZone* pZone_Ptr = pZone_.ptr();
        zs.pZone_ = pZone_Ptr->zones().append(pZone_Ptr);
    }
    else if (pZone_.valid())
    {
        zs.pZone_ = pZone_();
    }

    if (!cZone_.empty())
    {
        cellZone* cZone_Ptr = cZone_.ptr();
        zs.cZone_ = cZone_Ptr->zones().append(cZone_Ptr);
    }
    else if (cZone_.valid())
    {
        zs.cZone_ = cZone_();
    }

    if (!fZone_.empty())
    {
        faceZone* fZone_Ptr = fZone_.ptr();
        zs.fZone_ = fZone_Ptr->zones().append(fZone_Ptr);
    }
    else if (fZone_.valid())
    {
        zs.fZone_ = fZone_();
    }

    return zs;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::zoneSet::operator=(const zoneSet& zs)
{
    if (!zs.pZone_.empty())
    {
        pZone_ = zs.pZone_;
    }

    if (!zs.cZone_.empty())
    {
        cZone_ = zs.cZone_;
    }

    if (!zs.fZone_.empty())
    {
        fZone_ = zs.fZone_;
    }
}


// ************************************************************************* //
