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

#include "transformer.H"
#include "diagTensor.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline Foam::transformer Foam::transformer::translation(const vector& t)
{
    return transformer(t, true, tensor::I, false, false);
}


inline Foam::transformer Foam::transformer::scaling(const tensor& T)
{
    return transformer(vector::zero, false, T, true, false);
}


inline Foam::transformer Foam::transformer::rotation(const tensor& T)
{
    return transformer(vector::zero, false, T, false, true);
}


// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

inline Foam::transformer::transformer
(
    const vector& t,
    const bool translates,
    const tensor& T,
    const bool scales,
    const bool rotates
)
:
    t_(t),
    translates_(translates),
    T_(T),
    scales_(scales),
    rotates_(rotates)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::transformer::transformer()
:
    t_(Zero),
    translates_(false),
    T_(tensor::I),
    scales_(false),
    rotates_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::transformer::t() const
{
    return t_;
}


inline bool Foam::transformer::translates() const
{
    return translates_ && !transforms();
}


inline const Foam::tensor& Foam::transformer::T() const
{
    return T_;
}


inline Foam::tensor Foam::transformer::invT() const
{
    if (scales_ && rotates_)
    {
        return inv(T());
    }
    else if (scales_)
    {
        return inv(diag(T()));
    }
    else if (rotates_)
    {
        return T().T();
    }
    else
    {
        return tensor::I;
    }
}


inline bool Foam::transformer::scales() const
{
    return scales_ && !rotates_;
}


inline bool Foam::transformer::rotates() const
{
    return !scales_ && rotates_;
}


inline bool Foam::transformer::transforms() const
{
    return scales_ || rotates_;
}


template<typename Type>
inline bool Foam::transformer::transforms() const
{
    return pTraits<Type>::rank != 0 && transforms();
}


inline bool Foam::transformer::transformsPosition() const
{
    return translates_ || transforms();
}


inline Foam::vector Foam::transformer::transformPosition
(
    const vector& p
) const
{
    if (translates_ && !transforms())
    {
        return p + t();
    }
    else if (!translates_ && transforms())
    {
        return T() & p;
    }
    else if (translates_ && transforms())
    {
        return (T() & p) + t();
    }
    else
    {
        return p;
    }
}


inline Foam::vector Foam::transformer::invTransformPosition
(
    const vector& p
) const
{
    if (translates_ && !transforms())
    {
        return p - t();
    }
    else if (!translates_ && transforms())
    {
        return invT() & p;
    }
    else if (translates_ && transforms())
    {
        return invT() & (p - t());
    }
    else
    {
        return p;
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline Foam::transformer Foam::inv(const transformer& tr)
{
    return transformer
    (
        tr.invT() & (-tr.t()),
        tr.translates_,
        tr.invT(),
        tr.scales_,
        tr.rotates_
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const transformer& tr1, const transformer& tr2)
{
    return (tr1.t() == tr2.t() && tr1.T() == tr2.T());
}


inline bool Foam::operator!=(const transformer& tr1, const transformer& tr2)
{
    return !operator==(tr1, tr2);
}


inline Foam::transformer Foam::operator&
(
    const transformer& tr1,
    const transformer& tr2
)
{
    return transformer
    (
        (tr1.T() & tr2.t()) + tr1.t(),
        tr1.translates_ || tr2.translates_,
        tr1.T() & tr2.T(),
        tr1.scales_ || tr2.scales_,
        tr1.rotates_ || tr2.rotates_
    );
}


// ************************************************************************* //
