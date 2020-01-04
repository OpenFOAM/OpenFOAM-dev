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
    false
);

const Foam::transformer Foam::transformer::I
(
    Zero,
    false,
    tensor::I,
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

    buf << '(' << s.t() << ',' << s.R() << ')';

    return buf.str();
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

    is  >> tr.t_ >> tr.R_ >> tr.rotates_;

    // Read end of transformer
    is.readEnd("transformer");

    // Check state of Istream
    is.check("operator>>(Istream&, transformer&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const transformer& tr)
{
    os  << token::BEGIN_LIST
        << tr.t() << token::SPACE << tr.R() << token::SPACE << tr.rotates()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
