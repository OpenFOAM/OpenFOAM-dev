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

#include "verbatimString.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::verbatimString::verbatimString(Istream& is)
:
    string()
{
    is >> *this;
}


Foam::Istream& Foam::operator>>(Istream& is, verbatimString& vs)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isVerbatimString())
    {
        vs = t.verbatimStringToken();
    }
    else if (t.isString())
    {
        vs = t.stringToken();
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected verbatimString, found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, verbatimString&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const verbatimString& w)
{
    os.write(w);
    os.check("Ostream& operator<<(Ostream&, const verbatimString&)");
    return os;
}


// ************************************************************************* //
