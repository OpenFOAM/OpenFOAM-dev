/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "variable.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::variable::variable(Istream& is)
:
    word()
{
    is >> *this;
}


Foam::Istream& Foam::operator>>(Istream& is, variable& v)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isVariable())
    {
        v = t.variableToken();
    }
    else if (t.isWord())
    {
        v = t.wordToken();
    }
    else if (t.isString())
    {
        // Convert string to word stripping invalid characters
        v = t.stringToken();
        string::stripInvalid<variable>(v);

        // flag empty strings and bad chars as an error
        if (v.empty() || v.size() != t.stringToken().size())
        {
            is.setBad();
            FatalIOErrorInFunction(is)
                << "wrong token type - expected word, found "
                   "non-word characters "
                << t.info()
                << exit(FatalIOError);
            return is;
        }
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected word, found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, word&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const variable& w)
{
    os.write(w);
    os.check("Ostream& operator<<(Ostream&, const variable&)");
    return os;
}


// ************************************************************************* //
