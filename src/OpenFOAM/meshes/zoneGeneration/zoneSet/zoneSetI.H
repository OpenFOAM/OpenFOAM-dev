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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ZoneType>
inline void Foam::zoneSet::checkValid(const tmpNrc<ZoneType>& zone) const
{
    if (!zone.valid())
    {
        FatalErrorInFunction
            << "zoneSet does not contain a " << ZoneType::typeName << nl;

        if (pZone_.valid())
        {
            FatalError
                << "    zoneSet contains pointZone " << pZone_().name() << nl;
        }

        if (cZone_.valid())
        {
            FatalError
                << "    zoneSet contains cellZone " << cZone_().name() << nl;
        }

        if (fZone_.valid())
        {
            FatalError
                << "    zoneSet contains facesZone " << fZone_().name() << nl;
        }

        FatalError << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::zoneSet::zoneSet()
:
    pZone_(nullptr),
    cZone_(nullptr),
    fZone_(nullptr)
{}


inline Foam::zoneSet::zoneSet(const zoneSet& zs)
:
    pZone_(zs.pZone_),
    cZone_(zs.cZone_),
    fZone_(zs.fZone_)
{}


inline Foam::zoneSet::zoneSet(pointZone* pointZonePtr)
:
    pZone_(pointZonePtr),
    cZone_(nullptr),
    fZone_(nullptr)
{}


inline Foam::zoneSet::zoneSet(cellZone* cellZonePtr)
:
    pZone_(nullptr),
    cZone_(cellZonePtr),
    fZone_(nullptr)
{}


inline Foam::zoneSet::zoneSet(faceZone* faceZonePtr)
:
    pZone_(nullptr),
    cZone_(nullptr),
    fZone_(faceZonePtr)
{}


inline Foam::zoneSet::zoneSet
(
    pointZone* pointZonePtr,
    cellZone* cellZonePtr,
    faceZone* faceZonePtr
)
:
    pZone_(pointZonePtr),
    cZone_(cellZonePtr),
    fZone_(faceZonePtr)
{}


inline Foam::zoneSet::zoneSet(const pointZone& pointZoneRef)
:
    pZone_(pointZoneRef),
    cZone_(nullptr),
    fZone_(nullptr)
{}


inline Foam::zoneSet::zoneSet(const cellZone& cellZoneRef)
:
    pZone_(nullptr),
    cZone_(cellZoneRef),
    fZone_(nullptr)
{}


inline Foam::zoneSet::zoneSet(const faceZone& faceZoneRef)
:
    pZone_(nullptr),
    cZone_(nullptr),
    fZone_(faceZoneRef)
{}


inline Foam::zoneSet::zoneSet(const zoneSet& zs, bool allowTransfer)
:
    pZone_(zs.pZone_, allowTransfer),
    cZone_(zs.cZone_, allowTransfer),
    fZone_(zs.fZone_, allowTransfer)
{}


inline Foam::zoneSet Foam::zoneSet::clone(const word& name) const
{
    return zoneSet
    (
        pZone_.valid()
      ? pZone_().clone(name).ptr()
      : nullptr,

        cZone_.valid()
      ? cZone_().clone(name).ptr()
      : nullptr,

        fZone_.valid()
      ? fZone_().clone(name).ptr()
      : nullptr
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::zoneSet::pValid() const
{

    return pZone_.valid();
}


inline bool Foam::zoneSet::cValid() const
{

    return cZone_.valid();
}


inline bool Foam::zoneSet::fValid() const
{

    return fZone_.valid();
}


inline bool Foam::zoneSet::valid() const
{

    return pZone_.valid() || cZone_.valid() || fZone_.valid();
}


template<>
inline bool Foam::zoneSet::valid<Foam::pointZone>() const
{
    return pValid();
}


template<>
inline bool Foam::zoneSet::valid<Foam::cellZone>() const
{
    return cValid();
}


template<>
inline bool Foam::zoneSet::valid<Foam::faceZone>() const
{
    return fValid();
}


inline const Foam::pointZone& Foam::zoneSet::pZone() const
{
    checkValid(pZone_);
    return pZone_();
}


inline const Foam::cellZone& Foam::zoneSet::cZone() const
{
    checkValid(cZone_);
    return cZone_();
}


inline const Foam::faceZone& Foam::zoneSet::fZone() const
{
    checkValid(fZone_);
    return fZone_();
}


template<>
inline const Foam::pointZone& Foam::zoneSet::zone<Foam::pointZone>() const
{
    return pZone();
}


template<>
inline const Foam::cellZone& Foam::zoneSet::zone<Foam::cellZone>() const
{
    return cZone();
}


template<>
inline const Foam::faceZone& Foam::zoneSet::zone<Foam::faceZone>() const
{
    return fZone();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline void Foam::zoneSet::operator=(const pointZone& pZone)
{
    pZone_ = pZone;
}


inline void Foam::zoneSet::operator=(const cellZone& cZone)
{
    cZone_ = cZone;
}


inline void Foam::zoneSet::operator=(const faceZone& fZone)
{
    fZone_ = fZone;
}


// ************************************************************************* //
