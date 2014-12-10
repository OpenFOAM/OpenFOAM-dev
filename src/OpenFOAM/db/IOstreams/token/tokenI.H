/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Clear any allocated storage (word or string)
inline void token::clear()
{
    if (type_ == WORD)
    {
        delete wordTokenPtr_;
    }
    else if (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING)
    {
        delete stringTokenPtr_;
    }
    else if (type_ == COMPOUND)
    {
        if (compoundTokenPtr_->okToDelete())
        {
            delete compoundTokenPtr_;
        }
        else
        {
            compoundTokenPtr_->refCount::operator--();
        }
    }

    type_ = UNDEFINED;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
inline token::token()
:
    type_(UNDEFINED),
    lineNumber_(0)
{}

// Construct as copy
inline token::token(const token& t)
:
    type_(t.type_),
    lineNumber_(t.lineNumber_)
{
    switch (type_)
    {
        case token::UNDEFINED:
        break;

        case PUNCTUATION:
            punctuationToken_ = t.punctuationToken_;
        break;

        case WORD:
            wordTokenPtr_ = new word(*t.wordTokenPtr_);
        break;

        case STRING:
        case VARIABLE:
        case VERBATIMSTRING:
            stringTokenPtr_ = new string(*t.stringTokenPtr_);
        break;

        case LABEL:
            labelToken_ = t.labelToken_;
        break;

        case FLOAT_SCALAR:
            floatScalarToken_ = t.floatScalarToken_;
        break;

        case DOUBLE_SCALAR:
            doubleScalarToken_ = t.doubleScalarToken_;
        break;

        case COMPOUND:
            compoundTokenPtr_ = t.compoundTokenPtr_;
            compoundTokenPtr_->refCount::operator++();
        break;

        case token::ERROR:
        break;
    }
}

// Construct punctuation character token
inline token::token(punctuationToken p, label lineNumber)
:
    type_(PUNCTUATION),
    punctuationToken_(p),
    lineNumber_(lineNumber)
{}

// Construct word token
inline token::token(const word& w, label lineNumber)
:
    type_(WORD),
    wordTokenPtr_(new word(w)),
    lineNumber_(lineNumber)
{}

// Construct string token
inline token::token(const string& s, label lineNumber)
:
    type_(STRING),
    stringTokenPtr_(new string(s)),
    lineNumber_(lineNumber)
{}

// Construct label token
inline token::token(const label l, label lineNumber)
:
    type_(LABEL),
    labelToken_(l),
    lineNumber_(lineNumber)
{}

// Construct floatScalar token
inline token::token(const floatScalar s, label lineNumber)
:
    type_(FLOAT_SCALAR),
    floatScalarToken_(s),
    lineNumber_(lineNumber)
{}

// Construct doubleScalar token
inline token::token(const doubleScalar s, label lineNumber)
:
    type_(DOUBLE_SCALAR),
    doubleScalarToken_(s),
    lineNumber_(lineNumber)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// Delete token clearing the storage used by word or string
inline token::~token()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline token::tokenType token::type() const
{
    return type_;
}

inline token::tokenType& token::type()
{
    return type_;
}

inline bool token::good() const
{
    return (type_ != ERROR && type_ != UNDEFINED);
}

inline bool token::undefined() const
{
    return (type_ == UNDEFINED);
}

inline bool token::error() const
{
    return (type_ == ERROR);
}

inline bool token::isPunctuation() const
{
    return (type_ == PUNCTUATION);
}

inline token::punctuationToken token::pToken() const
{
    if (type_ == PUNCTUATION)
    {
        return punctuationToken_;
    }
    else
    {
        parseError("punctuation character");
        return NULL_TOKEN;
    }
}

inline bool token::isWord() const
{
    return (type_ == WORD);
}

inline const word& token::wordToken() const
{
    if (type_ == WORD)
    {
        return *wordTokenPtr_;
    }
    else
    {
        parseError("word");
        return word::null;
    }
}

inline bool token::isVariable() const
{
    return (type_ == VARIABLE);
}

inline bool token::isString() const
{
    return (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING);
}

inline const string& token::stringToken() const
{
    if (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING)
    {
        return *stringTokenPtr_;
    }
    else
    {
        parseError("string");
        return string::null;
    }
}

inline bool token::isLabel() const
{
    return (type_ == LABEL);
}

inline label token::labelToken() const
{
    if (type_ == LABEL)
    {
        return labelToken_;
    }
    else
    {
        parseError("label");
        return 0;
    }
}

inline bool token::isFloatScalar() const
{
    return (type_ == FLOAT_SCALAR);
}

inline floatScalar token::floatScalarToken() const
{
    if (type_ == FLOAT_SCALAR)
    {
        return floatScalarToken_;
    }
    else
    {
        parseError("floatScalar");
        return 0.0;
    }
}


inline bool token::isDoubleScalar() const
{
    return (type_ == DOUBLE_SCALAR);
}

inline doubleScalar token::doubleScalarToken() const
{
    if (type_ == DOUBLE_SCALAR)
    {
        return doubleScalarToken_;
    }
    else
    {
        parseError("doubleScalar");
        return 0.0;
    }
}


inline bool token::isScalar() const
{
    return (type_ == FLOAT_SCALAR || type_ == DOUBLE_SCALAR);
}

inline scalar token::scalarToken() const
{
    if (type_ == FLOAT_SCALAR)
    {
        return floatScalarToken_;
    }
    else if (type_ == DOUBLE_SCALAR)
    {
        return doubleScalarToken_;
    }
    else
    {
        parseError("scalar");
        return 0.0;
    }
}

inline bool token::isNumber() const
{
    return (type_ == LABEL || isScalar());
}

inline scalar token::number() const
{
    if (type_ == LABEL)
    {
        return labelToken_;
    }
    else if (isScalar())
    {
        return scalarToken();
    }
    else
    {
        parseError("number (label or scalar)");
        return 0.0;
    }
}

inline bool token::isCompound() const
{
    return (type_ == COMPOUND);
}

inline const token::compound& token::compoundToken() const
{
    if (type_ == COMPOUND)
    {
        return *compoundTokenPtr_;
    }
    else
    {
        parseError("compound");
        return *compoundTokenPtr_;
    }
}


inline label token::lineNumber() const
{
    return lineNumber_;
}

inline label& token::lineNumber()
{
    return lineNumber_;
}


inline void token::setBad()
{
    clear();
    type_ = ERROR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void token::operator=(const token& t)
{
    clear();
    type_ = t.type_;

    switch (type_)
    {
        case token::UNDEFINED:
        break;

        case PUNCTUATION:
            punctuationToken_ = t.punctuationToken_;
        break;

        case WORD:
            wordTokenPtr_ = new word(*t.wordTokenPtr_);
        break;

        case STRING:
        case VARIABLE:
        case VERBATIMSTRING:
            stringTokenPtr_ = new string(*t.stringTokenPtr_);
        break;

        case LABEL:
            labelToken_ = t.labelToken_;
        break;

        case FLOAT_SCALAR:
            floatScalarToken_ = t.floatScalarToken_;
        break;

        case DOUBLE_SCALAR:
            doubleScalarToken_ = t.doubleScalarToken_;
        break;

        case COMPOUND:
            compoundTokenPtr_ = t.compoundTokenPtr_;
            compoundTokenPtr_->refCount::operator++();
        break;

        case token::ERROR:
        break;
    }

    lineNumber_ = t.lineNumber_;
}

inline void token::operator=(const punctuationToken p)
{
    clear();
    type_ = PUNCTUATION;
    punctuationToken_ = p;
}

inline void token::operator=(word* wPtr)
{
    clear();
    type_ = WORD;
    wordTokenPtr_ = wPtr;
}

inline void token::operator=(const word& w)
{
    operator=(new word(w));
}

inline void token::operator=(string* sPtr)
{
    clear();
    type_ = STRING;
    stringTokenPtr_ = sPtr;
}

inline void token::operator=(const string& s)
{
    operator=(new string(s));
}

inline void token::operator=(const label l)
{
    clear();
    type_ = LABEL;
    labelToken_ = l;
}

inline void token::operator=(const floatScalar s)
{
    clear();
    type_ = FLOAT_SCALAR;
    floatScalarToken_ = s;
}

inline void token::operator=(const doubleScalar s)
{
    clear();
    type_ = DOUBLE_SCALAR;
    doubleScalarToken_ = s;
}

inline void token::operator=(token::compound* cPtr)
{
    clear();
    type_ = COMPOUND;
    compoundTokenPtr_ = cPtr;
}


inline bool token::operator==(const token& t) const
{
    if (type_ != t.type_)
    {
        return false;
    }

    switch (type_)
    {
        case token::UNDEFINED:
            return true;

        case PUNCTUATION:
            return punctuationToken_ == t.punctuationToken_;

        case WORD:
            return *wordTokenPtr_ == *t.wordTokenPtr_;

        case STRING:
        case VARIABLE:
        case VERBATIMSTRING:
            return *stringTokenPtr_ == *t.stringTokenPtr_;

        case LABEL:
            return labelToken_ == t.labelToken_;

        case FLOAT_SCALAR:
            return equal(floatScalarToken_, t.floatScalarToken_);

        case DOUBLE_SCALAR:
            return equal(doubleScalarToken_, t.doubleScalarToken_);

        case COMPOUND:
            return compoundTokenPtr_ == t.compoundTokenPtr_;

        case token::ERROR:
            return true;
    }

    return false;
}

inline bool token::operator==(const punctuationToken p) const
{
    return (type_ == PUNCTUATION && punctuationToken_ == p);
}

inline bool token::operator==(const word& w) const
{
    return (type_ == WORD && wordToken() == w);
}

inline bool token::operator==(const string& s) const
{
    return
    (
        (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING)
     && stringToken() == s
    );
}

inline bool token::operator==(const label l) const
{
    return (type_ == LABEL && labelToken_ == l);
}

inline bool token::operator==(const floatScalar s) const
{
    return (type_ == FLOAT_SCALAR && equal(floatScalarToken_, s));
}

inline bool token::operator==(const doubleScalar s) const
{
    return (type_ == DOUBLE_SCALAR && equal(doubleScalarToken_, s));
}

inline bool token::operator!=(const token& t) const
{
    return !operator==(t);
}

inline bool token::operator!=(const punctuationToken p) const
{
    return !operator==(p);
}

inline bool token::operator!=(const word& w) const
{
    return !operator==(w);
}

inline bool token::operator!=(const string& s) const
{
    return !operator==(s);
}

inline bool token::operator!=(const floatScalar s) const
{
    return !operator==(s);
}

inline bool token::operator!=(const doubleScalar s) const
{
    return !operator==(s);
}

inline bool token::operator!=(const label l) const
{
    return !operator==(l);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
