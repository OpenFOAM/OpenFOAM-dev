/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "WallCollisionRecord.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::scalar Foam::WallCollisionRecord<Type>::errorCosAngle(1.0 + 1e-6);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::WallCollisionRecord<Type>::WallCollisionRecord()
:
    accessed_(false),
    pRel_(),
    data_(Zero)
{}


template<class Type>
Foam::WallCollisionRecord<Type>::WallCollisionRecord
(
    bool accessed,
    const vector& pRel,
    const Type& data
)
:
    accessed_(accessed),
    pRel_(pRel),
    data_(data)
{}


template<class Type>
Foam::WallCollisionRecord<Type>::WallCollisionRecord
(
    const WallCollisionRecord<Type>& wCR
)
:
    accessed_(wCR.accessed_),
    pRel_(wCR.pRel_),
    data_(wCR.data_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::WallCollisionRecord<Type>::~WallCollisionRecord()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::WallCollisionRecord<Type>::operator=
(
    const WallCollisionRecord<Type>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    accessed_ = rhs.accessed_;
    pRel_ = rhs.pRel_;
    data_ = rhs.data_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "WallCollisionRecordIO.C"


// ************************************************************************* //
