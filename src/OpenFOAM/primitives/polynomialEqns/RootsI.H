/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

template <Foam::direction N>
inline Foam::Roots<N>::Roots()
:
    types_(0)
{
    forAll(*this, i)
    {
        type(i, roots::nan);
    }
}


template <Foam::direction N>
inline Foam::Roots<N>::Roots(const roots::type t, const scalar x)
:
    types_(0)
{
    forAll(*this, i)
    {
        this->v_[i] = x;
        type(i, t);
    }
}


template <Foam::direction N>
inline Foam::Roots<N>::Roots
(
    const roots::type t,
    const scalar x,
    const Roots<N - 1>& xs
)
:
    types_(0)
{
    this->v_[0] = x;
    type(0, t);
    forAll(xs, i)
    {
        this->v_[i+1] = xs[i];
        type(i + 1, xs.type(i));
    }
}


template <Foam::direction N>
inline Foam::Roots<N>::Roots
(
    const Roots<N - 1>& xs,
    const roots::type t,
    const scalar x
)
:
    types_(0)
{
    forAll(xs, i)
    {
        this->v_[i] = xs[i];
        type(i, xs.type(i));
    }
    this->v_[N-1] = x;
    type(N - 1, t);
}


template <Foam::direction N>
template <Foam::direction M>
inline Foam::Roots<N>::Roots
(
    const Roots<M>& xs,
    const Roots<N - M>& ys
)
:
    types_(0)
{
    forAll(xs, i)
    {
        this->v_[i] = xs[i];
        type(i, xs.type(i));
    }
    forAll(ys, i)
    {
        this->v_[i + M] = ys[i];
        type(i + M, ys.type(i));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <Foam::direction N>
inline void Foam::Roots<N>::type
(
    const direction i,
    const roots::type t
)
{
    types_ += (t - type(i)) << 3*i;
}


template <Foam::direction N>
inline Foam::roots::type Foam::Roots<N>::type(const direction i) const
{
    return static_cast<roots::type>((types_ >> 3*i) & 7);
}


// ************************************************************************* //
