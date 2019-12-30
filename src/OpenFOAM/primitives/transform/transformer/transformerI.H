/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    R_(sphericalTensor::I),
    hasR_(false)
{}


inline Foam::transformer::transformer
(
    const vector& t,
    const tensor& R,
    bool hasR
)
:
    t_(t),
    R_(R),
    hasR_(hasR)
{}


inline Foam::transformer::transformer(const vector& t)
:
    t_(t),
    R_(sphericalTensor::I),
    hasR_(false)
{}


inline Foam::transformer::transformer(const tensor& R)
:
    t_(Zero),
    R_(R),
    hasR_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::transformer::t() const
{
    return t_;
}


inline const Foam::tensor& Foam::transformer::R() const
{
    return R_;
}


inline bool Foam::transformer::hasR() const
{
    return hasR_;
}


inline Foam::vector& Foam::transformer::t()
{
    return t_;
}


inline Foam::tensor& Foam::transformer::R()
{
    // Assume that non-const access to R changes it from I, so set
    // hasR to true

    hasR_ = true;

    return R_;
}


inline Foam::vector Foam::transformer::transformPosition
(
    const vector& v
) const
{
    if (hasR_)
    {
        return t() + (R() & v);
    }
    else
    {
        return t() + v;
    }
}


inline Foam::pointField Foam::transformer::transformPosition
(
    const pointField& pts
) const
{
    tmp<pointField> tfld;

    if (hasR_)
    {
        tfld = t() + (R() & pts);
    }
    else
    {
        tfld = t() + pts;
    }
    return tfld();
}


inline Foam::vector Foam::transformer::invTransformPosition
(
    const vector& v
) const
{
    if (hasR_)
    {
        return (R().T() & (v - t()));
    }
    else
    {
        return v - t();
    }
}


inline Foam::pointField Foam::transformer::invTransformPosition
(
    const pointField& pts
) const
{
    tmp<pointField> tfld;

    if (hasR_)
    {
        tfld = (R().T() & (pts - t()));
    }
    else
    {
        tfld = pts - t();
    }
    return tfld();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::transformer::operator&=
(
    const transformer& tr
)
{
    t_ += tr.t_;
    R_ = tr.R_ & R_;

    // If either of the two objects has hasR_ as true, then inherit
    // it, otherwise, these should both be I tensors.
    hasR_ = tr.hasR_ || hasR_;
}


inline void Foam::transformer::operator=(const vector& t)
{
    t_ = t;
}


inline void Foam::transformer::operator+=(const vector& t)
{
    t_ += t;
}


inline void Foam::transformer::operator-=(const vector& t)
{
    t_ -= t;
}


inline void Foam::transformer::operator=(const tensor& R)
{
    hasR_ = true;

    R_ = R;
}


inline void Foam::transformer::operator&=(const tensor& R)
{
    hasR_ = true;

    R_ = R & R_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline Foam::transformer Foam::inv(const transformer& tr)
{
    return transformer(-tr.t(), tr.R().T(), tr.hasR());
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==
(
    const transformer& tr1,
    const transformer& tr2
)
{
    return (tr1.t() == tr2.t() && tr1.R() == tr2.R());
}


inline bool Foam::operator!=
(
    const transformer& tr1,
    const transformer& tr2
)
{
    return !operator==(tr1, tr2);
}


inline Foam::transformer Foam::operator+
(
    const transformer& tr,
    const vector& t
)
{
    return transformer(tr.t() + t, tr.R(), tr.hasR());
}


inline Foam::transformer Foam::operator+
(
    const vector& t,
    const transformer& tr
)
{
    return transformer(t + tr.t(), tr.R(), tr.hasR());
}


inline Foam::transformer Foam::operator-
(
    const transformer& tr,
    const vector& t
)
{
    return transformer(tr.t() - t, tr.R(), tr.hasR());
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
        tr1.R() & tr2.R(),
        (tr1.hasR() || tr2.hasR())
    );
}


// ************************************************************************* //
