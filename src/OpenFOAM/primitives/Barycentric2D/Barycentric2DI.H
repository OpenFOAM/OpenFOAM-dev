/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Barycentric2D<Cmpt>::Barycentric2D()
{}


template<class Cmpt>
inline Foam::Barycentric2D<Cmpt>::Barycentric2D(const Foam::zero)
:
    Barycentric2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::Barycentric2D<Cmpt>::Barycentric2D
(
    const Cmpt& va,
    const Cmpt& vb,
    const Cmpt& vc
)
{
    this->v_[A] = va;
    this->v_[B] = vb;
    this->v_[C] = vc;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::Barycentric2D<Cmpt>::a() const
{
    return this->v_[A];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric2D<Cmpt>::b() const
{
    return this->v_[B];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric2D<Cmpt>::c() const
{
    return this->v_[C];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric2D<Cmpt>::a()
{
    return this->v_[A];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric2D<Cmpt>::b()
{
    return this->v_[B];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric2D<Cmpt>::c()
{
    return this->v_[C];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt operator&
(
    const Barycentric2D<Cmpt>& b1,
    const Barycentric2D<Cmpt>& b2
)
{
    return b1.a()*b2.a() + b1.b()*b2.b() + b1.c()*b2.c();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
