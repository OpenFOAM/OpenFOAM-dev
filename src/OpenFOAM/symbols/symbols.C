/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "symbols.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::symbols::tokeniser::push(const token& t)
{
    label end = (start_+size_)%tokens_.size();
    tokens_[end] = t;
    if (size_ == tokens_.size())
    {
        start_ = tokens_.fcIndex(start_);
    }
    else
    {
        size_++;
    }
}


Foam::token Foam::symbols::tokeniser::pop()
{
    token t = tokens_[start_];
    start_ = tokens_.fcIndex(start_);
    --size_;
    return t;
}


void Foam::symbols::tokeniser::unpop(const token& t)
{
    ++size_;
    start_ = tokens_.rcIndex(start_);
    tokens_[start_] = t;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::symbols::tokeniser::tokeniser(Istream& is)
:
    is_(is),
    tokens_(100),
    start_(0),
    size_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::symbols::tokeniser::stream()
{
    return is_;
}


bool Foam::symbols::tokeniser::hasToken() const
{
    return size_ || is_.good();
}


Foam::token Foam::symbols::tokeniser::nextToken()
{
    if (size_ == 0)
    {
        token t(is_);

        if (t.isWord())
        {
            splitWord(t.wordToken());
            return pop();
        }
        else
        {
            return t;
        }
    }
    else
    {
        return pop();
    }
}


void Foam::symbols::tokeniser::putBack(const token& t)
{
    if (size_ == 0)
    {
        push(t);
    }
    else
    {
        unpop(t);
    }
}


void Foam::symbols::tokeniser::splitWord(const word& w)
{
    size_t start = 0;
    for (size_t i=0; i<w.size(); ++i)
    {
        if (!valid(w[i]))
        {
            if (i > start)
            {
                word subWord = w(start, i-start);
                if (isdigit(subWord[0]) || subWord[0] == token::SUBTRACT)
                {
                    push(token(readScalar(IStringStream(subWord)())));
                }
                else
                {
                    push(token(subWord));
                }
            }
            if (w[i] != token::SPACE)
            {
                if (isdigit(w[i]))
                {
                    push(token(readScalar(IStringStream(w[i])())));
                }
                else
                {
                    push(token::punctuationToken(w[i]));
                }
            }
            start = i+1;
        }
    }
    if (start < w.size())
    {
        word subWord = w(start, w.size()-start);
        if (isdigit(subWord[0]) || subWord[0] == token::SUBTRACT)
        {
            push(token(readScalar(IStringStream(subWord)())));
        }
        else
        {
            push(token(subWord));
        }
    }
}


bool Foam::symbols::tokeniser::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // divide
     && c != ';'   // end statement
     && c != '{'   // begin sub-dictionary
     && c != '}'   // end sub-dictionary
     && c != '('   // begin expression
     && c != ')'   // end expression
     && c != '['   // begin dimensions/units
     && c != ']'   // end dimensions/units
     && c != ':'   // separate dimensions/units
     && c != '^'   // power
     && c != '*'   // multiply
    );
}


Foam::label Foam::symbols::tokeniser::priority(const token& t)
{
    if (!t.isPunctuation())
    {
        return 0;
    }
    else if
    (
        t.pToken() == token::MULTIPLY
     || t.pToken() == token::DIVIDE
    )
    {
        return 2;
    }
    else if (t.pToken() == '^')
    {
        return 3;
    }
    else
    {
        return 0;
    }
}


// ************************************************************************* //
