/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Ostream::decrIndent()
{
    if (indentLevel_ == 0)
    {
        cerr<< "Ostream::decrIndent() : attempt to decrement 0 indent level"
            << std::endl;
    }
    else
    {
        indentLevel_--;
    }
}


Foam::Ostream& Foam::Ostream::writeKeyword(const keyType& kw)
{
    return Foam::writeKeyword(*this, kw);
}


Foam::Ostream& Foam::Ostream::write(const token& t)
{
    switch (t.type())
    {
        case token::UNDEFINED:
            write("UNDEFINED");
            WarningInFunction << "Undefined token" << Foam::endl;
        break;

        case token::PUNCTUATION:
            // Cast to char required to work around a bug in gcc versions < 10
            write(static_cast<char>(t.pToken()));
            // write(t.pToken());
        break;

        case token::WORD:
            write(t.wordToken());
        break;

        case token::FUNCTIONNAME:
            write(t.functionNameToken());
        break;

        case token::VARIABLE:
            write(t.variableToken());
        break;

        case token::STRING:
            write(t.stringToken());
        break;

        case token::VERBATIMSTRING:
            write(t.verbatimStringToken());
        break;

        case token::INTEGER_32:
            write(t.integer32Token());
        break;

        case token::INTEGER_64:
            write(t.integer64Token());
        break;

        case token::UNSIGNED_INTEGER_32:
            write(t.unsignedInteger32Token());
        break;

        case token::UNSIGNED_INTEGER_64:
            write(t.unsignedInteger64Token());
        break;

        case token::FLOAT_SCALAR:
            write(t.floatScalarToken());
        break;

        case token::DOUBLE_SCALAR:
            write(t.doubleScalarToken());
        break;

        case token::LONG_DOUBLE_SCALAR:
            write(t.longDoubleScalarToken());
        break;

        case token::COMPOUND:
            *this << t.compoundToken();
        break;

        case token::ERROR:
            write("ERROR");
            WarningInFunction << "Error token" << Foam::endl;
        break;

        default:
            write("UNKNOWN");
            SeriousErrorInFunction << "Unknown token" << Foam::endl;
    }

    // Check state of stream
    check("Ostream& Ostream::write(const token&)");

    return *this;
}


Foam::Ostream& Foam::Ostream::writeCompoundTag(const word& typeName)
{
    if (token::compound::isCompound(typeName))
    {
        write(typeName);
        write(token::SPACE);
    }

    return *this;
}


// ************************************************************************* //
