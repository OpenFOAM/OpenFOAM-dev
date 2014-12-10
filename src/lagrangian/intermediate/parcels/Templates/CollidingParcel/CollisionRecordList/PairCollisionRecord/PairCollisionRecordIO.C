/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord(Istream& is)
:
    origProcOfOther_(readLabel(is)),
    origIdOfOther_(readLabel(is)),
    data_(is)
{
    // Check state of Istream
    is.check
    (
        "Foam::PairCollisionRecord<Type>::PairCollisionRecord(Foam::Istream&)"
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, PairCollisionRecord<Type>& pCR)
{
    is  >> pCR.origProcOfOther_ >> pCR.origIdOfOther_ >> pCR.data_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream&"
        "Foam::operator>>(Foam::Istream&, Foam::PairCollisionRecord<Type>&)"
    );

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PairCollisionRecord<Type>& pCR
)
{
    os  << pCR.origProcOfOther_
        << token::SPACE << pCR.origIdOfOther_
        << token::SPACE << pCR.data_;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::PairCollisionRecord<Type>&)"
    );

    return os;
}


// ************************************************************************* //
