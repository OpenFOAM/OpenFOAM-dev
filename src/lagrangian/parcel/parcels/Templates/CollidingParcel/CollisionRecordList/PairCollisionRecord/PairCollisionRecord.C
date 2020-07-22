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

#include "PairCollisionRecord.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord()
:
    origProcOfOther_(0),
    origIdOfOther_(-1),
    data_(Zero)
{}


template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord
(
    bool accessed,
    label origProcOfOther,
    label origIdOfOther,
    const Type& data
)
:
    origProcOfOther_(origProcOfOther + 1),
    origIdOfOther_(origIdOfOther),
    data_(data)
{
    // Default assignment to origProcOfOther_ assumes accessed is true

    if (!accessed)
    {
        setUnaccessed();
    }
}


template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord
(
    const PairCollisionRecord<Type>& pCR
)
:
    origProcOfOther_(pCR.origProcOfOther_),
    origIdOfOther_(pCR.origIdOfOther_),
    data_(pCR.data_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PairCollisionRecord<Type>::~PairCollisionRecord()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::PairCollisionRecord<Type>::operator=
(
    const PairCollisionRecord<Type>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    origProcOfOther_ = rhs.origProcOfOther_;
    origIdOfOther_ = rhs.origIdOfOther_;
    data_ = rhs.data_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "PairCollisionRecordIO.C"


// ************************************************************************* //
