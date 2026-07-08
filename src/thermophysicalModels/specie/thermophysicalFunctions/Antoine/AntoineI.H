/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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

#include "Antoine.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Function1Type>
inline Foam::Function1s::AntoineFunction1<Function1Type>::AntoineFunction1
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, Function1Type>(name),
    AntoineCoeffs(Function1Type::typeName, units, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Function1Type>
inline Foam::Function1s::AntoineFunction1<Function1Type>::~AntoineFunction1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::Function1s::Antoine::value(const scalar T) const
{
    return pow(scalar(10), A_ - B_/(C_ + T));
}


inline Foam::scalar Foam::Function1s::logAntoine::value(const scalar T) const
{
    return ln10_*(A_ - B_/(C_ + T));
}


inline Foam::scalar Foam::Function1s::inverseAntoine::value
(
    const scalar p
) const
{
    return B_/(A_ - log10(p)) - C_;
}


inline Foam::scalar Foam::Function1s::Antoine::derivative(const scalar T) const
{
    return ln10_*value(T)*B_/sqr(C_ + T);
}


inline Foam::scalar Foam::Function1s::logAntoine::derivative
(
    const scalar T
) const
{
    return ln10_*B_/sqr(C_ + T);
}


inline Foam::scalar Foam::Function1s::inverseAntoine::derivative
(
    const scalar p
) const
{
    return B_/(ln10_*p*sqr(A_ - log10(p)));
}


inline Foam::scalar Foam::Function1s::Antoine::integral
(
    const scalar T1,
    const scalar T2
) const
{
    NotImplemented;
    return NaN;
}


inline Foam::scalar Foam::Function1s::logAntoine::integral
(
    const scalar T1,
    const scalar T2
) const
{
    NotImplemented;
    return NaN;
}


inline Foam::scalar Foam::Function1s::inverseAntoine::integral
(
    const scalar p1,
    const scalar p2
) const
{
    NotImplemented;
    return NaN;
}


template<class Function1Type>
inline void Foam::Function1s::AntoineFunction1<Function1Type>::write
(
    Ostream& os,
    const unitSets& units
) const
{
    AntoineCoeffs::write(os, units);
}


// ************************************************************************* //
