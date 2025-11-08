/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "error.H"
#include "token.H"

#include "IOstreams.H"
#include "OTstream.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::token::token(Istream& is)
:
    type_(UNDEFINED)
{
    is.read(*this);
}


// * * * * * * * * * * * * IOstream operators  * * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, token& t)
{
    t.clear();
    return is.read(t);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token& t)
{
    return os.write(t);
}


ostream& Foam::operator<<(ostream& os, const token::punctuationToken& pt)
{
    return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::punctuationToken& pt)
{
    return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::compound& ct)
{
    os << ct.type() << token::SPACE;
    ct.write(os);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ostream& Foam::operator<<(ostream& os, const InfoProxy<token>& ip)
{
    const token& t = ip.t_;

    os  << "on line " << t.lineNumber();

    switch (t.type())
    {
        case token::UNDEFINED:
            os  << " an undefined token";
        break;

        case token::PUNCTUATION:
            os  << " the punctuation token " << '\'' << t.pToken() << '\'';
        break;

        case token::WORD:
            os  << " the word " << '\'' << t.wordToken() << '\'';
        break;

        case token::STRING:
            os  << " the string " << t.stringToken();
        break;

        case token::VERBATIMSTRING:
            os  << " the verbatim string " << t.verbatimStringToken();
        break;

        case token::FUNCTIONNAME:
            os  << " the functionName " << t.functionNameToken();
        break;

        case token::VARIABLE:
            os  << " the variable " << t.variableToken();
        break;

        case token::INTEGER_32:
            os  << " the 32-bit integer " << t.integer32Token();
        break;

        case token::INTEGER_64:
            os  << " the 64-bit integer " << t.integer64Token();
        break;

        case token::UNSIGNED_INTEGER_32:
            os  << " the unsigned 32-bit integer "
                << t.unsignedInteger32Token();
        break;

        case token::UNSIGNED_INTEGER_64:
            os  << " the unsigned 64-bit integer "
                << t.unsignedInteger64Token();
        break;

        case token::FLOAT_SCALAR:
            os  << " the floatScalar " << t.floatScalarToken();
        break;

        case token::DOUBLE_SCALAR:
            os  << " the doubleScalar " << t.doubleScalarToken();
        break;

        case token::LONG_DOUBLE_SCALAR:
            os  << " the longDoubleScalar " << t.longDoubleScalarToken();
        break;

        case token::COMPOUND:
        {
            if (t.compoundToken().empty())
            {
                os  << " the empty compound of type "
                    << t.compoundToken().type();
            }
            else
            {
                os  << " the compound of type "
                    << t.compoundToken().type();
            }
        }
        break;

        case token::ERROR:
            os  << " an error";
        break;

        default:
            os  << " an unknown token type " << '\'' << int(t.type()) << '\'';
    }

    return os;
}


template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<token>& ip)
{
    const token& t = ip.t_;

    os  << "on line " << t.lineNumber();

    switch (t.type())
    {
        case token::UNDEFINED:
            os  << " an undefined token";
        break;

        case token::PUNCTUATION:
            os  << " the punctuation token " << '\'' << t.pToken() << '\'';
        break;

        case token::WORD:
            os  << " the word " << '\'' << t.wordToken() << '\'';
        break;

        case token::STRING:
            os  << " the string " << t.stringToken();
        break;

        case token::VERBATIMSTRING:
            os  << " the verbatim string " << t.verbatimStringToken();
        break;

        case token::FUNCTIONNAME:
            os  << " the functionName " << t.functionNameToken();
        break;

        case token::VARIABLE:
            os  << " the variable " << t.variableToken();
        break;

        break;
        case token::INTEGER_32:
            os  << " the 32-bit integer " << t.integer32Token();
        break;

        case token::INTEGER_64:
            os  << " the 64-bit integer " << t.integer64Token();
        break;

        case token::UNSIGNED_INTEGER_32:
            os  << " the unsigned 32-bit integer "
                << t.unsignedInteger32Token();
        break;

        case token::UNSIGNED_INTEGER_64:
            os  << " the unsigned 64-bit integer "
                << t.unsignedInteger64Token();
        break;

        case token::FLOAT_SCALAR:
            os  << " the floatScalar " << t.floatScalarToken();
        break;

        case token::DOUBLE_SCALAR:
            os  << " the doubleScalar " << t.doubleScalarToken();
        break;

        case token::LONG_DOUBLE_SCALAR:
            os  << " the longDoubleScalar " << t.longDoubleScalarToken();
        break;

        case token::COMPOUND:
        {
            if (t.compoundToken().empty())
            {
                os  << " the empty compound of type "
                    << t.compoundToken().type();
            }
            else
            {
                os  << " the compound of type "
                    << t.compoundToken().type();
            }
        }
        break;

        case token::ERROR:
            os  << " an error";
        break;

        default:
            os  << " an unknown token type "  << '\'' << int(t.type()) << '\'';
    }

    return os;
}


// ************************************************************************* //
