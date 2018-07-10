/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
inline Foam::Barycentric<Cmpt>::Barycentric()
{}


template<class Cmpt>
inline Foam::Barycentric<Cmpt>::Barycentric(const Foam::zero)
:
    Barycentric::vsType(Zero)
{}


template<class Cmpt>
inline Foam::Barycentric<Cmpt>::Barycentric
(
    const Cmpt& va,
    const Cmpt& vb,
    const Cmpt& vc,
    const Cmpt& vd
)
{
    this->v_[A] = va;
    this->v_[B] = vb;
    this->v_[C] = vc;
    this->v_[D] = vd;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::a() const
{
    return this->v_[A];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::b() const
{
    return this->v_[B];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::c() const
{
    return this->v_[C];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::d() const
{
    return this->v_[D];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::a()
{
    return this->v_[A];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::b()
{
    return this->v_[B];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::c()
{
    return this->v_[C];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::d()
{
    return this->v_[D];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt operator&(const Barycentric<Cmpt>& b1, const Barycentric<Cmpt>& b2)
{
    return b1.a()*b2.a() + b1.b()*b2.b() + b1.c()*b2.c() + b1.d()*b2.d();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
