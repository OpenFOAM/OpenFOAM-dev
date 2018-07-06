/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "cellShape.H"
#include "token.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, cellShape& s)
{
    bool readEndBracket = false;

    // Read the 'name' token for the symbol
    token t(is);

    if (t.isPunctuation())
    {
        if (t.pToken() == token::BEGIN_LIST)
        {
            readEndBracket = true;

            is >> t;
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << t.info()
                << exit(FatalIOError);
        }
    }

    // it is allowed to have either a word or a number describing the model
    if (t.isLabel())
    {
        s.m = cellModeller::lookup(int(t.labelToken()));
    }
    else if (t.isWord())
    {
        s.m = cellModeller::lookup(t.wordToken());
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Bad type of token for cellShape symbol " << t.info()
            << exit(FatalIOError);
        return is;
    }

    // Check that a model was found
    if (!s.m)
    {
        FatalIOErrorInFunction(is)
            << "CellShape has unknown model " << t.info()
            << exit(FatalIOError);
        return is;
    }

    // Read the geometry labels
    is >> static_cast<labelList&>(s);

    if (readEndBracket)
    {
        // Read end)
        is.readEnd("cellShape");
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const cellShape & s)
{
    // Write beginning of record
    os << token::BEGIN_LIST;

    // Write the list label for the symbol (ONE OR THE OTHER !!!)
    os << (s.m)->index() << token::SPACE;

    // Write the model name instead of the label (ONE OR THE OTHER !!!)
    // os << (s.m)->name() << token::SPACE;

    // Write the geometry
    os << static_cast<const labelList&>(s);

    // End of record
    os << token::END_LIST;

    return os;
}


template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<cellShape>& ip)
{
    const cellShape& cs = ip.t_;

    if (isNull(cs.model()))
    {
        os  << "    cellShape has no model!\n";
    }
    else
    {
        os  << cs.model().info() << endl;
    }

    os  << "\tGeom:\tpoint\tlabel\txyz\n";

    forAll(cs, i)
    {
        os  << "\t\t" << i << "\t" << cs[i] << endl;
    }

    return os;
}


// ************************************************************************* //
