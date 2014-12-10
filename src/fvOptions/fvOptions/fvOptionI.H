/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::fv::option::name() const
{
    return name_;
}


inline const Foam::fvMesh& Foam::fv::option::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::fv::option::coeffs() const
{
    return coeffs_;
}


inline bool Foam::fv::option::active() const
{
    return active_;
}


inline Foam::scalar Foam::fv::option::timeStart() const
{
    return timeStart_;
}


inline Foam::scalar Foam::fv::option::duration() const
{
    return duration_;
}


inline bool Foam::fv::option::inTimeLimits(const scalar time) const
{
    return
    (
        (timeStart_ < 0)
     ||
        (
            (mesh_.time().value() >= timeStart_)
         && (mesh_.time().value() <= (timeStart_ + duration_))
        )
    );
}


inline const Foam::fv::option::selectionModeType&
Foam::fv::option::selectionMode() const
{
    return selectionMode_;
}


inline const Foam::word& Foam::fv::option::cellSetName() const
{
    return cellSetName_;
}


inline Foam::scalar Foam::fv::option::V() const
{
    return V_;
}


inline const Foam::labelList& Foam::fv::option::cells() const
{
    return cells_;
}


inline void Foam::fv::option::setApplied(const label fieldI)
{
    applied_[fieldI] = true;
}


inline bool& Foam::fv::option::active()
{
    return active_;
}


inline Foam::scalar& Foam::fv::option::timeStart()
{
    return timeStart_;
}


inline Foam::scalar& Foam::fv::option::duration()
{
    return duration_;
}


inline const Foam::word& Foam::fv::option::nbrRegionName() const
{
    return nbrRegionName_;
}


inline const Foam::meshToMesh& Foam::fv::option::meshInterp() const
{
    if (!meshInterpPtr_.valid())
    {
        FatalErrorIn("const meshToMesh& meshInterp() const")
            << "Interpolation object not set"
            << abort(FatalError);
    }

    return meshInterpPtr_();
}


// ************************************************************************* //
