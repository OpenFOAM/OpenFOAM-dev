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

Description
    Istream constructor and IOstream operators for keyType.

\*---------------------------------------------------------------------------*/

#include "keyType.H"
#include "regExp.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::keyType Foam::keyType::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::keyType::keyType(Istream& is)
:
    word(),
    isPattern_(false)
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
    if (literalMatch || !isPattern_)
    {
        // check as string
        return (str == *this);
    }
    else
    {
        // check as regex
        return regExp(*this).match(str);
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

    if (t.isWord())
    {
        kw = t.wordToken();
    }
    else if (t.isString())
    {
        // Assign from string. Set as regular expression.
        kw = t.stringToken();
        kw.isPattern_ = true;

        // flag empty strings as an error
        if (kw.empty())
        {
            is.setBad();
            FatalIOErrorIn("operator>>(Istream&, keyType&)", is)
                << "empty word/expression "
                << exit(FatalIOError);
            return is;
        }
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, keyType&)", is)
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
    os.write(kw);
    os.check("Ostream& operator<<(Ostream&, const keyType&)");
    return os;
}


// ************************************************************************* //
