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

#include "scalarRange.H"
#include "token.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::scalarRange::debug(::Foam::debug::debugSwitch("scalarRange", 0));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarRange::scalarRange()
:
    type_(EMPTY),
    value_(0),
    value2_(0)
{}


Foam::scalarRange::scalarRange(const scalar lower, const scalar upper)
:
    type_(RANGE),
    value_(lower),
    value2_(upper)
{
    // mark invalid range as empty
    if (lower > upper)
    {
        type_ = EMPTY;
        value_ = value2_ = 0;
    }
}


Foam::scalarRange::scalarRange(Istream& is)
:
    type_(EXACT),
    value_(0),
    value2_(0)
{
    is >> *this;

    if (scalarRange::debug)
    {
        Info<<"constructed scalarRange: " << *this << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::scalarRange::empty() const
{
    return type_ == EMPTY;
}


bool Foam::scalarRange::valid() const
{
    return type_ != EMPTY;
}


bool Foam::scalarRange::isExact() const
{
    return type_ == EXACT;
}


Foam::scalar Foam::scalarRange::value() const
{
    return value_;
}


Foam::scalar Foam::scalarRange::lower() const
{
    if (type_ == UPPER)
    {
        return -Foam::GREAT;
    }
    else
    {
        return value_;
    }
}

Foam::scalar Foam::scalarRange::upper() const
{
    switch (type_)
    {
        case LOWER:
            return Foam::GREAT;
            break;

        case RANGE:
            return value2_;
            break;

        default:
            return value_;
            break;
    }
}


bool Foam::scalarRange::selected(const scalar value) const
{
    switch (type_)
    {
        case LOWER:
            return value >= value_;

        case UPPER:
            return value <= value_;

        case RANGE:
            return value >= value_ && value <= value2_;

        case EXACT:
            return value == value_;

        default:
            return false;
    }
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::scalarRange::operator==(const scalarRange& range) const
{
    return
    (
        type_ == range.type_
     && value_ == range.value_
     && value2_ == range.value2_
    );
}


bool Foam::scalarRange::operator!=(const scalarRange& range) const
{
    return !(operator==(range));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, scalarRange& range)
{
    range.type_ = scalarRange::EXACT;
    range.value_ = 0;
    range.value2_ = 0;

    List<token> toks(4);
    label nTok = 0;

    // skip leading ','
    do
    {
        is.read(toks[nTok]);
        is.check("scalarRange token read");
    }
    while
    (
        toks[nTok].isPunctuation()
     && toks[nTok].pToken() == token::COMMA
    );

    ++nTok;

    // looks like ':VALUE'
    if
    (
        toks[nTok-1].isPunctuation()
     && toks[nTok-1].pToken() == token::COLON
    )
    {
        range.type_ = scalarRange::UPPER;
        is.read(toks[nTok++]);
        is.check("scalarRange token read");
    }

    // a number is now required
    if (!toks[nTok-1].isNumber())
    {
        is.setBad();
        range.type_ = scalarRange::EMPTY;
        range.value_ = range.value2_ = 0;
        Info<< "rejected ill-formed or empty range:";
        for (label i=0; i<nTok; ++i)
        {
            Info<< " " << toks[i];
        }
        Info<< endl;
        return is;
    }

    range.value_ = toks[nTok-1].number();
    is.read(toks[nTok++]);
    is.check("scalarRange token read");

    if (scalarRange::debug)
    {
        Info<<"tokens:";
        for (label i=0; i<nTok; ++i)
        {
            Info<< " " << toks[i];
        }
        Info<< endl;
    }

    // could be 'VALUE:' or 'VALUE:VALUE'
    if
    (
        toks[nTok-1].isPunctuation()
     && toks[nTok-1].pToken() == token::COLON
    )
    {
        if (range.type_ == scalarRange::UPPER)
        {
            is.setBad();
            range.type_ = scalarRange::EMPTY;
            range.value_ = range.value2_ = 0;
            Info<< "rejected ill-formed range:";
            for (label i=0; i<nTok; ++i)
            {
                Info<< " " << toks[i];
            }
            Info<< endl;
            return is;
        }

        is.read(toks[nTok++]);
        is.check("scalarRange token read");

        if (scalarRange::debug)
        {
            Info<<"tokens:";
            for (label i=0; i<nTok; ++i)
            {
                Info<< " " << toks[i];
            }
            Info<< endl;
        }


        // if there is a number, we have 'VALUE:VALUE' and not simply 'VALUE:'
        if (toks[nTok-1].isNumber())
        {
            range.type_ = scalarRange::RANGE;
            range.value2_ = toks[nTok-1].number();
            is.read(toks[nTok++]);
            is.check("scalarRange token read");
        }
        else
        {
            range.type_ = scalarRange::LOWER;
        }
    }

    if (scalarRange::debug)
    {
        Info<<"tokens:";
        for (label i=0; i<nTok; ++i)
        {
            Info<< " " << toks[i];
        }
        Info<< endl;
    }


    // some remaining tokens, but they are not the next comma
    // - this is a problem!
    if
    (
        toks[nTok-1].good()
     &&
        (
            !toks[nTok-1].isPunctuation()
         || toks[nTok-1].pToken() != token::COMMA
        )
    )
    {
        is.setBad();
        range.type_ = scalarRange::EMPTY;
        range.value_ = range.value2_ = 0;

        Info<< "rejected ill-formed range:";
        for (label i=0; i<nTok; ++i)
        {
            Info<< " " << toks[i];
        }
        Info<< endl;
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const scalarRange& range)
{
    switch (range.type_)
    {
        case scalarRange::LOWER:
            os << range.value_ << " <=> Inf";
            break;

        case scalarRange::UPPER:
            os << "-Inf <=> " << range.value_;
            break;

        case scalarRange::RANGE:
            os << range.value_ << " <=> " << range.value2_;
            break;

        case scalarRange::EXACT:
            os << range.value_;
            break;

        default:
            os << "empty";
            break;
    }

    return os;
}


// ************************************************************************* //
