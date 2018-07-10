/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::WallCollisionRecord<Type>::match
(
    const vector& pRel,
    scalar radius
)
{
    scalar magpRel_= mag(pRel_);

    scalar magpRel = mag(pRel);

    // Using the new data as the acceptance criterion
    scalar cosAcceptanceAngle = magpRel/radius;

    if (cosAcceptanceAngle > errorCosAngle)
    {
        Info<< "pRel_ " << pRel_ << " " << magpRel_ << nl
            << "pRel " << pRel << " " << magpRel << nl
            << "unit vector dot product " << (pRel & pRel_)/(magpRel_*magpRel)
            << nl << "cosAcceptanceAngle " << cosAcceptanceAngle
            << endl;

        FatalErrorInFunction
            << "Problem with matching WallCollisionRecord." << nl
            << "The given radius, " << radius << ", is smaller than distance "
            << "to the relative position of the WallInteractionSite, "
            << magpRel << nl
            << abort(FatalError);
    }

    // Are the test and recorded pRel (relative position vectors)
    // aligned to within the calculated tolerance?
    bool matched = (pRel & pRel_)/(magpRel_*magpRel) > cosAcceptanceAngle;

    if (matched)
    {
        pRel_ = pRel;
    }

    return matched;
}


template<class Type>
inline const Foam::vector&
Foam::WallCollisionRecord<Type>::pRel() const
{
    return pRel_;
}


template<class Type>
inline const Type&
Foam::WallCollisionRecord<Type>::collisionData() const
{
    return data_;
}


template<class Type>
inline Type& Foam::WallCollisionRecord<Type>::collisionData()
{
    return data_;
}


template<class Type>
inline bool Foam::WallCollisionRecord<Type>::accessed() const
{
    return accessed_;
}


template<class Type>
inline void Foam::WallCollisionRecord<Type>::setAccessed()
{
    accessed_ = true;
}


template<class Type>
inline void Foam::WallCollisionRecord<Type>::setUnaccessed()
{
    accessed_ = false;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::operator==
(
    const WallCollisionRecord<Type>& a,
    const WallCollisionRecord<Type>& b
)
{
    return
    (
        a.accessed_ == b.accessed_
     && a.pRel_ == b.pRel_
     && a.data_ == b.data_
    );
}


template<class Type>
inline bool Foam::operator!=
(
    const WallCollisionRecord<Type>& a,
    const WallCollisionRecord<Type>& b
)
{
    return !(a == b);
}


// ************************************************************************* //
