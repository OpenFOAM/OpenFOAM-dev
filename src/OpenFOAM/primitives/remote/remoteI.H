/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "remote.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::remote::remote()
:
    proci(-1),
    elementi(-1)
{}


inline Foam::remote::remote(const label p, const label e)
:
    proci(p),
    elementi(e)
{}


inline Foam::remote::remote(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const remote& a, const remote& b)
{
    return a.proci == b.proci && a.elementi == b.elementi;
}


inline bool Foam::operator!=(const remote& a, const remote& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const remote& p)
{
    return os << p.proci << token::SPACE << p.elementi;
}


inline Foam::Istream& Foam::operator>>(Istream& is, remote& p)
{
    return is >> p.proci >> p.elementi;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

template<class Cmpt>
class Tensor;

inline remote transform(const Tensor<scalar>&, const remote& p)
{
    return p;
}

}


// ************************************************************************* //
