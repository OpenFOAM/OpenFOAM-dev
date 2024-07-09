/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "tokenTupleN.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tokenTupleN::tokenTupleN()
:
    tokens_(),
    offsets_()
{}


Foam::tokenTupleN::tokenTupleN(Istream& is)
:
    tokens_(4),
    offsets_(2)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tokenTupleN::~tokenTupleN()
{}


// * * * * * * * * * * * * * * IOStream Operators * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, tokenTupleN& ttn)
{
    is.readBegin("StringTupleN");

    label level = 0;

    while (is.good())
    {
        if (level == 0)
        {
            ttn.offsets_.append(ttn.tokens_.size());
        }

        const token t(is);

        if (t == token::BEGIN_LIST)
        {
            level ++;
        }

        if (t == token::END_LIST)
        {
            if (level != 0)
            {
                level --;
            }
            else
            {
                is.putBack(t);
                is.readEnd("StringTupleN");
                break;
            }
        }

        ttn.tokens_.append(t);
    }

    // Check state of Istream
    is.check("operator>>(Istream& is, tokenTupleN&)");

    return is;
}


// ************************************************************************* //
