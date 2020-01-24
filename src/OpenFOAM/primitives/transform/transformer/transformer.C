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
#include "IOstreams.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::transformer::typeName =
    "transformer";

const Foam::transformer Foam::transformer::zero
(
    Zero,
    false,
    Zero,
    false,
    false
);

const Foam::transformer Foam::transformer::I
(
    Zero,
    false,
    tensor::I,
    false,
    false
);

const Foam::transformer Foam::transformer::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transformer::transformer(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::name(const transformer& s)
{
    OStringStream buf;

    buf << '(' << s.t() << ',' << s.T() << ')';

    return buf.str();
}


void Foam::transformer::transformPosition
(
    pointField& res,
    const pointField& pts
) const
{
    if (translates_ && !transforms())
    {
        res = pts + t();
    }
    else if (!translates_ && transforms())
    {
        res = T() & pts;
    }
    else if (translates_ && transforms())
    {
        res = (T() & pts) + t();
    }
}


Foam::tmp<Foam::pointField> Foam::transformer::transformPosition
(
    const pointField& pts
) const
{
    if (translates_ && !transforms())
    {
        return pts + t();
    }
    else if (!translates_ && transforms())
    {
        return T() & pts;
    }
    else if (translates_ && transforms())
    {
        return (T() & pts) + t();
    }
    else
    {
        return pts;
    }
}


Foam::tmp<Foam::pointField> Foam::transformer::invTransformPosition
(
    const pointField& pts
) const
{
    if (translates_ && !transforms())
    {
        return pts - t();
    }
    else if (!translates_ && transforms())
    {
        return T().T() & pts;
    }
    else if (translates_ && transforms())
    {
        return (T().T() & (pts - t()));
    }
    else
    {
        return pts;
    }
}


void Foam::transformer::invTransformPosition
(
    pointField& res,
    const pointField& pts
) const
{
    if (translates_ && !transforms())
    {
        res = pts - t();
    }
    else if (!translates_ && transforms())
    {
        res = T().T() & pts;
    }
    else if (translates_ && transforms())
    {
        res = (T().T() & (pts - t()));
    }
}


template<>
Foam::tmp<Foam::Field<bool>> Foam::transformer::transform
(
    const Field<bool>& fld
) const
{
    return fld;
}


template<>
Foam::tmp<Foam::Field<bool>> Foam::transformer::transform
(
    const tmp<Field<bool>>& tfld
) const
{
    return tfld;
}


template<>
Foam::tmp<Foam::Field<Foam::label>> Foam::transformer::transform
(
    const Field<label>& fld
) const
{
    return fld;
}


template<>
Foam::tmp<Foam::Field<Foam::label>> Foam::transformer::transform
(
    const tmp<Field<label>>& tfld
) const
{
    return tfld;
}


template<>
Foam::tmp<Foam::Field<Foam::scalar>> Foam::transformer::transform
(
    const Field<scalar>& fld
) const
{
    return fld;
}


template<>
Foam::tmp<Foam::Field<Foam::scalar>> Foam::transformer::transform
(
    const tmp<Field<scalar>>& tfld
) const
{
    return tfld;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, transformer& tr)
{
    // Read beginning of transformer
    is.readBegin("transformer");

    is  >> tr.translates_ >> tr.t_ >> tr.scales_ >> tr.rotates_ >> tr.T_;

    // Read end of transformer
    is.readEnd("transformer");

    // Check state of Istream
    is.check("operator>>(Istream&, transformer&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const transformer& tr)
{
    os  << token::BEGIN_LIST
        << tr.translates_ << token::SPACE << tr.t_ << token::SPACE
        << tr.scales_ << token::SPACE << tr.rotates_ << token::SPACE << tr.T_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
