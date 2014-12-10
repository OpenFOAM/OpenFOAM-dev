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

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline specie::specie
(
    const word& name,
    const scalar nMoles,
    const scalar molWeight
)
:
    name_(name),
    nMoles_(nMoles),
    molWeight_(molWeight)
{}


inline specie::specie
(
    const scalar nMoles,
    const scalar molWeight
)
:
    nMoles_(nMoles),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline specie::specie(const specie& st)
:
    name_(st.name_),
    nMoles_(st.nMoles_),
    molWeight_(st.molWeight_)
{}


inline specie::specie(const word& name, const specie& st)
:
    name_(name),
    nMoles_(st.nMoles_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& specie::name() const
{
    return name_;
}


inline scalar specie::W() const
{
    return molWeight_;
}


inline scalar specie::nMoles() const
{
    return nMoles_;
}


inline scalar specie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void specie::operator=(const specie& st)
{
    //name_ = st.name_;
    nMoles_ = st.nMoles_;
    molWeight_ = st.molWeight_;
}


inline void specie::operator+=(const specie& st)
{
    scalar sumNmoles = max(nMoles_ + st.nMoles_, SMALL);

    molWeight_ =
        nMoles_/sumNmoles*molWeight_
      + st.nMoles_/sumNmoles*st.molWeight_;

    nMoles_ = sumNmoles;
}


inline void specie::operator-=(const specie& st)
{
    scalar diffnMoles = nMoles_ - st.nMoles_;
    if (mag(diffnMoles) < SMALL)
    {
        diffnMoles = SMALL;
    }

    molWeight_ =
        nMoles_/diffnMoles*molWeight_
      - st.nMoles_/diffnMoles*st.molWeight_;

    nMoles_ = diffnMoles;
}


inline void specie::operator*=(const scalar s)
{
    nMoles_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline specie operator+(const specie& st1, const specie& st2)
{
    scalar sumNmoles = max(st1.nMoles_ + st2.nMoles_, SMALL);

    return specie
    (
        sumNmoles,
        st1.nMoles_/sumNmoles*st1.molWeight_
      + st2.nMoles_/sumNmoles*st2.molWeight_
    );
}


inline specie operator-(const specie& st1, const specie& st2)
{
    scalar diffNmoles = st1.nMoles_ - st2.nMoles_;
    if (mag(diffNmoles) < SMALL)
    {
        diffNmoles = SMALL;
    }

    return specie
    (
        diffNmoles,
        st1.nMoles_/diffNmoles*st1.molWeight_
      - st2.nMoles_/diffNmoles*st2.molWeight_
    );
}


inline specie operator*(const scalar s, const specie& st)
{
    return specie
    (
        s*st.nMoles_,
        st.molWeight_
    );
}


inline specie operator==(const specie& st1, const specie& st2)
{
    return st2 - st1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
