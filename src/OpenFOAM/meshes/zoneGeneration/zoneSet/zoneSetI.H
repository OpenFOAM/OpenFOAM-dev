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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::zoneSet::zoneSet()
:
    pZone(nullptr),
    cZone(nullptr),
    fZone(nullptr)
{}


inline Foam::zoneSet::zoneSet(const zoneSet& zs)
:
    pZone(zs.pZone),
    cZone(zs.cZone),
    fZone(zs.fZone)
{}


// inline Foam::zoneSet::zoneSet(zoneSet&& zs)
// :
//     pZone(std::move(zs.pZone)),
//     cZone(std::move(zs.cZone)),
//     fZone(std::move(zs.fZone))
// {}


// inline Foam::zoneSet::zoneSet(const zoneSet&& zs)
// :
//     pZone(std::move(zs.pZone)),
//     cZone(std::move(zs.cZone)),
//     fZone(std::move(zs.fZone))
// {}


inline Foam::zoneSet::zoneSet(pointZone* pointZonePtr)
:
    pZone(pointZonePtr),
    cZone(nullptr),
    fZone(nullptr)
{}


inline Foam::zoneSet::zoneSet(cellZone* cellZonePtr)
:
    pZone(nullptr),
    cZone(cellZonePtr),
    fZone(nullptr)
{}


inline Foam::zoneSet::zoneSet(faceZone* faceZonePtr)
:
    pZone(nullptr),
    cZone(nullptr),
    fZone(faceZonePtr)
{}


inline Foam::zoneSet::zoneSet
(
    pointZone* pointZonePtr,
    cellZone* cellZonePtr,
    faceZone* faceZonePtr
)
:
    pZone(pointZonePtr),
    cZone(cellZonePtr),
    fZone(faceZonePtr)
{}


inline Foam::zoneSet::zoneSet(const pointZone& pointZoneRef)
:
    pZone(pointZoneRef),
    cZone(nullptr),
    fZone(nullptr)
{}


inline Foam::zoneSet::zoneSet(const cellZone& cellZoneRef)
:
    pZone(nullptr),
    cZone(cellZoneRef),
    fZone(nullptr)
{}


inline Foam::zoneSet::zoneSet(const faceZone& faceZoneRef)
:
    pZone(nullptr),
    cZone(nullptr),
    fZone(faceZoneRef)
{}


inline Foam::zoneSet::zoneSet(const zoneSet& zs, bool allowTransfer)
:
    pZone(zs.pZone, allowTransfer),
    cZone(zs.cZone, allowTransfer),
    fZone(zs.fZone, allowTransfer)
{}


inline Foam::zoneSet Foam::zoneSet::clone(const word& name) const
{
    return zoneSet
    (
        pZone.valid()
      ? pZone().clone(name).ptr()
      : nullptr,

        cZone.valid()
      ? cZone().clone(name).ptr()
      : nullptr,

        fZone.valid()
      ? fZone().clone(name).ptr()
      : nullptr
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::zoneSet::valid() const
{

    return pZone.valid() || cZone.valid() || fZone.valid();
}


// ************************************************************************* //
