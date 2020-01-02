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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::transformer::transformer()
:
    t_(Zero),
    translates_(false),
    R_(tensor::I),
    rotates_(false)
{}


inline Foam::transformer::transformer(const vector& t)
:
    t_(t),
    translates_(true),
    R_(tensor::I),
    rotates_(false)
{}


inline Foam::transformer::transformer(const tensor& R)
:
    t_(Zero),
    translates_(false),
    R_(R),
    rotates_(true)
{}


inline Foam::transformer::transformer(const vector& t, const tensor& R)
:
    t_(t),
    translates_(true),
    R_(R),
    rotates_(true)
{}


inline Foam::transformer::transformer
(
    const vector& t,
    const bool translates,
    const tensor& R,
    const bool rotates
)
:
    t_(t),
    translates_(translates),
    R_(R),
    rotates_(rotates)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::transformer::t() const
{
    return t_;
}


inline bool Foam::transformer::translates() const
{
    return translates_;
}


inline const Foam::tensor& Foam::transformer::R() const
{
    return R_;
}


inline bool Foam::transformer::rotates() const
{
    return rotates_;
}


inline bool Foam::transformer::transforms() const
{
    return rotates_;
}


template<typename Type>
inline bool Foam::transformer::transforms() const
{
    return pTraits<Type>::rank != 0 && (translates_ || rotates_);
}


inline bool Foam::transformer::transformsPosition() const
{
    return translates_ || rotates_;
}


inline Foam::vector& Foam::transformer::t()
{
    // Assume that non-const access to t changes it from zero
    translates_ = true;

    return t_;
}


inline Foam::tensor& Foam::transformer::R()
{
    // Assume that non-const access to R changes it from I
    rotates_ = true;

    return R_;
}


inline Foam::vector Foam::transformer::transformPosition
(
    const vector& v
) const
{
    if (translates_ && !rotates_)
    {
        return v + t();
    }
    else if (!translates_ && rotates_)
    {
        return R() & v;
    }
    else if (translates_ && rotates_)
    {
        return (R() & v) + t();
    }
    else
    {
        return v;
    }
}


inline Foam::tmp<Foam::pointField> Foam::transformer::transformPosition
(
    const pointField& pts
) const
{
    if (translates_ && !rotates_)
    {
        return pts + t();
    }
    else if (!translates_ && rotates_)
    {
        return R() & pts;
    }
    else if (translates_ && rotates_)
    {
        return (R() & pts) + t();
    }
    else
    {
        return pts;
    }
}


inline Foam::vector Foam::transformer::invTransformPosition
(
    const vector& v
) const
{
    if (translates_ && !rotates_)
    {
        return v - t();
    }
    else if (!translates_ && rotates_)
    {
        return R().T() & v;
    }
    else if (translates_ && rotates_)
    {
        return (R().T() & (v - t()));
    }
    else
    {
        return v;
    }
}


inline Foam::tmp<Foam::pointField> Foam::transformer::invTransformPosition
(
    const pointField& pts
) const
{
    if (translates_ && !rotates_)
    {
        return pts - t();
    }
    else if (!translates_ && rotates_)
    {
        return R().T() & pts;
    }
    else if (translates_ && rotates_)
    {
        return (R().T() & (pts - t()));
    }
    else
    {
        return pts;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::transformer::operator&=
(
    const transformer& tr
)
{
    t_ += tr.t_;

    // If either of the two objects has translates_ as true, then inherit
    // it, otherwise, these should both be zero vectors.
    translates_ = tr.translates_ || translates_;

    R_ = tr.R_ & R_;

    // If either of the two objects has rotates_ as true, then inherit
    // it, otherwise, these should both be I tensors.
    rotates_ = tr.rotates_ || rotates_;
}


inline void Foam::transformer::operator=(const vector& t)
{
    translates_ = true;
    t_ = t;
}


inline void Foam::transformer::operator+=(const vector& t)
{
    translates_ = true;
    t_ += t;
}


inline void Foam::transformer::operator-=(const vector& t)
{
    translates_ = true;
    t_ -= t;
}


inline void Foam::transformer::operator=(const tensor& R)
{
    rotates_ = true;
    R_ = R;
}


inline void Foam::transformer::operator&=(const tensor& R)
{
    rotates_ = true;
    R_ = R & R_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline Foam::transformer Foam::inv(const transformer& tr)
{
    if (tr.translates() && !tr.rotates())
    {
        return transformer(-tr.t());
    }
    else if (!tr.translates() && tr.rotates())
    {
        return transformer(tr.R().T());
    }
    else if (tr.translates() && tr.rotates())
    {
        return transformer(tr.R().T() & (-tr.t()), tr.R().T());
    }
    else
    {
        return transformer();
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const transformer& tr1, const transformer& tr2)
{
    return (tr1.t() == tr2.t() && tr1.R() == tr2.R());
}


inline bool Foam::operator!=(const transformer& tr1, const transformer& tr2)
{
    return !operator==(tr1, tr2);
}


inline Foam::transformer Foam::operator+(const transformer& tr, const vector& t)
{
    return transformer(tr.t() + t, true, tr.R(), tr.rotates());
}


inline Foam::transformer Foam::operator+(const vector& t, const transformer& tr)
{
    return transformer(t + tr.t(), true, tr.R(), tr.rotates());
}


inline Foam::transformer Foam::operator-(const transformer& tr, const vector& t)
{
    return transformer(tr.t() - t, true, tr.R(), tr.rotates());
}


inline Foam::transformer Foam::operator&
(
    const transformer& tr1,
    const transformer& tr2
)
{
    return transformer
    (
        tr1.t() + tr2.t(),
        tr1.translates() || tr2.translates(),
        tr1.R() & tr2.R(),
        tr1.rotates() || tr2.rotates()
    );
}


// ************************************************************************* //
