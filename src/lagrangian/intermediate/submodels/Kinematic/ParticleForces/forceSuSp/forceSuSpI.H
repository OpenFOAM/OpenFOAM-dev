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

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

inline Foam::forceSuSp::forceSuSp()
{}


inline Foam::forceSuSp::forceSuSp
(
    const Tuple2<vector, scalar>& fs
)
:
    Tuple2<vector, scalar>(fs)
{}


inline Foam::forceSuSp::forceSuSp(const vector& Su, const scalar Sp)
{
    first() = Su;
    second() = Sp;
}


inline Foam::forceSuSp::forceSuSp(Istream& is)
:
    Tuple2<vector, scalar>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::forceSuSp::Su() const
{
    return first();
}


inline Foam::scalar Foam::forceSuSp::Sp() const
{
    return second();
}


inline Foam::vector& Foam::forceSuSp::Su()
{
    return first();
}


inline Foam::scalar& Foam::forceSuSp::Sp()
{
    return second();
}


// * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * * //

inline void Foam::forceSuSp::operator=(const forceSuSp& susp)
{
    first() = susp.first();
    second() = susp.second();
}


inline void Foam::forceSuSp::operator+=(const forceSuSp& susp)
{
    first() += susp.first();
    second() += susp.second();
}


inline void Foam::forceSuSp::operator-=(const forceSuSp& susp)
{
    first() -= susp.first();
    second() -= susp.second();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::forceSuSp Foam::operator+
(
    const forceSuSp& susp1,
    const forceSuSp& susp2
)
{
    return forceSuSp
    (
        susp1.first() + susp2.first(),
        susp1.second() + susp2.second()
    );
}


inline Foam::forceSuSp Foam::operator*
(
    const scalar s,
    const forceSuSp& susp
)
{
    return forceSuSp(susp.first()*s, susp.second()*s);
}


// ************************************************************************* //
