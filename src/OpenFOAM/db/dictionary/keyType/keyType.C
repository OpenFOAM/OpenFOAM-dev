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

#include "keyType.H"
#include "regExp.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::keyType Foam::keyType::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::keyType::keyType(const token& t)
:
    variable(),
    type_(UNDEFINED)
{
    operator=(t);
}


Foam::keyType::keyType(Istream& is)
:
    variable(),
    type_(UNDEFINED)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::keyType::match
(
    const std::string& str,
    bool literalMatch
) const
{
    if (literalMatch || !isPattern())
    {
        // Check as string
        return (str == *this);
    }
    else
    {
        // Check as regex
        return regExp(*this).match(str);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::keyType::operator=(const token& t)
{
    if (t.isWord())
    {
        operator=(t.wordToken());
    }
    else if (t.isFunctionName())
    {
        operator=(t.functionNameToken());
    }
    else if (t.isVariable())
    {
        operator=(t.variableToken());
    }
    else if (t.isString())
    {
        // Assign from string. Set as pattern.
        operator=(t.stringToken());

        // An empty pattern string is a fatal error
        if (empty())
        {
            FatalErrorInFunction
                << "Empty pattern string"
                << exit(FatalIOError);
        }
    }
    else
    {
        variable::clear();
        type_ = UNDEFINED;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, keyType& kw)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    kw = t;

    if (kw.isUndefined())
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected word or string, found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, keyType&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const keyType& kw)
{
    os.writeQuoted(kw, kw.isPattern());
    os.check("Ostream& operator<<(Ostream&, const keyType&)");
    return os;
}


// * * * * * * * * * * * * * * * * Global Function * * * * * * * * * * * * * //

Foam::Ostream& Foam::writeKeyword(Foam::Ostream& os, const keyType& kw)
{
    //- Indentation of the entry from the start of the keyword
    static const unsigned short entryIndentation_ = 16;

    os.indent();
    os << kw;

    label nSpaces = entryIndentation_ - label(kw.size());

    // Pattern is surrounded by quotes
    if (kw.isPattern())
    {
        nSpaces -= 2;
    }

    nSpaces = max(nSpaces, 1);

    while (nSpaces--)
    {
        os.write(char(token::SPACE));
    }

    return os;
}


// ************************************************************************* //
