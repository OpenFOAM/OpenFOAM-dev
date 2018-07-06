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

#include <iostream>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::token::clear()
{
    if (type_ == WORD)
    {
        delete wordTokenPtr_;
    }
    else if (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING)
    {
        delete stringTokenPtr_;
    }
    else if (type_ == LONG_DOUBLE_SCALAR)
    {
        delete longDoubleScalarTokenPtr_;
    }
    else if (type_ == COMPOUND)
    {
        if (compoundTokenPtr_->unique())
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

inline Foam::token::token()
:
    type_(UNDEFINED),
    lineNumber_(0)
{}


inline Foam::token::token(const token& t)
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

        case LONG_DOUBLE_SCALAR:
            longDoubleScalarTokenPtr_ =
                new longDoubleScalar(*t.longDoubleScalarTokenPtr_);
        break;

        case COMPOUND:
            compoundTokenPtr_ = t.compoundTokenPtr_;
            compoundTokenPtr_->refCount::operator++();
        break;

        case token::ERROR:
        break;
    }
}


inline Foam::token::token(punctuationToken p, label lineNumber)
:
    type_(PUNCTUATION),
    punctuationToken_(p),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const word& w, label lineNumber)
:
    type_(WORD),
    wordTokenPtr_(new word(w)),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const string& s, label lineNumber)
:
    type_(STRING),
    stringTokenPtr_(new string(s)),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const label l, label lineNumber)
:
    type_(LABEL),
    labelToken_(l),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const floatScalar s, label lineNumber)
:
    type_(FLOAT_SCALAR),
    floatScalarToken_(s),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const doubleScalar s, label lineNumber)
:
    type_(DOUBLE_SCALAR),
    doubleScalarToken_(s),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const longDoubleScalar s, label lineNumber)
:
    type_(LONG_DOUBLE_SCALAR),
    longDoubleScalarTokenPtr_(new longDoubleScalar(s)),
    lineNumber_(lineNumber)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::token::~token()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::token::tokenType  Foam::token::type() const
{
    return type_;
}

inline Foam::token::tokenType&  Foam::token::type()
{
    return type_;
}

inline bool Foam::token::good() const
{
    return (type_ != ERROR && type_ != UNDEFINED);
}

inline bool Foam::token::undefined() const
{
    return (type_ == UNDEFINED);
}

inline bool Foam::token::error() const
{
    return (type_ == ERROR);
}

inline bool Foam::token::isPunctuation() const
{
    return (type_ == PUNCTUATION);
}

inline Foam::token::punctuationToken  Foam::token::pToken() const
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

inline bool Foam::token::isWord() const
{
    return (type_ == WORD);
}

inline const Foam::word& Foam::token::wordToken() const
{
    if (type_ == WORD)
    {
        return *wordTokenPtr_;
    }
    else
    {
        parseError(word::typeName);
        return word::null;
    }
}

inline bool Foam::token::isVariable() const
{
    return (type_ == VARIABLE);
}

inline bool Foam::token::isString() const
{
    return (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING);
}

inline const Foam::string& Foam::token::stringToken() const
{
    if (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING)
    {
        return *stringTokenPtr_;
    }
    else
    {
        parseError(string::typeName);
        return string::null;
    }
}

inline bool Foam::token::isLabel() const
{
    return (type_ == LABEL);
}

inline Foam::label Foam::token::labelToken() const
{
    if (type_ == LABEL)
    {
        return labelToken_;
    }
    else
    {
        parseError(pTraits<label>::typeName);
        return 0;
    }
}

inline bool Foam::token::isFloatScalar() const
{
    return (type_ == FLOAT_SCALAR);
}

inline Foam::floatScalar Foam::token::floatScalarToken() const
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


inline bool Foam::token::isDoubleScalar() const
{
    return (type_ == DOUBLE_SCALAR);
}

inline Foam::doubleScalar Foam::token::doubleScalarToken() const
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


inline bool Foam::token::isLongDoubleScalar() const
{
    return (type_ == LONG_DOUBLE_SCALAR);
}

inline Foam::longDoubleScalar Foam::token::longDoubleScalarToken() const
{
    if (type_ == LONG_DOUBLE_SCALAR)
    {
        return *longDoubleScalarTokenPtr_;
    }
    else
    {
        parseError("longDoubleScalar");
        return 0.0;
    }
}


inline bool Foam::token::isScalar() const
{
    return
    (
        type_ == FLOAT_SCALAR
     || type_ == DOUBLE_SCALAR
     || type_ == LONG_DOUBLE_SCALAR
    );
}

inline Foam::scalar Foam::token::scalarToken() const
{
    if (type_ == FLOAT_SCALAR)
    {
        return floatScalarToken_;
    }
    else if (type_ == DOUBLE_SCALAR)
    {
        return doubleScalarToken_;
    }
    else if (type_ == LONG_DOUBLE_SCALAR)
    {
        return *longDoubleScalarTokenPtr_;
    }
    else
    {
        parseError(pTraits<scalar>::typeName);
        return 0.0;
    }
}

inline bool Foam::token::isNumber() const
{
    return (type_ == LABEL || isScalar());
}

inline Foam::scalar Foam::token::number() const
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

inline bool Foam::token::isCompound() const
{
    return (type_ == COMPOUND);
}

inline const Foam::token::compound& Foam::token::compoundToken() const
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


inline Foam::label Foam::token::lineNumber() const
{
    return lineNumber_;
}

inline Foam::label& Foam::token::lineNumber()
{
    return lineNumber_;
}


inline void Foam::token::setBad()
{
    clear();
    type_ = ERROR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::token::operator=(const token& t)
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

        case LONG_DOUBLE_SCALAR:
            longDoubleScalarTokenPtr_ =
                new longDoubleScalar(*t.longDoubleScalarTokenPtr_);
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

inline void Foam::token::operator=(const punctuationToken p)
{
    clear();
    type_ = PUNCTUATION;
    punctuationToken_ = p;
}

inline void Foam::token::operator=(word* wPtr)
{
    clear();
    type_ = WORD;
    wordTokenPtr_ = wPtr;
}

inline void Foam::token::operator=(const word& w)
{
    operator=(new word(w));
}

inline void Foam::token::operator=(string* sPtr)
{
    clear();
    type_ = STRING;
    stringTokenPtr_ = sPtr;
}

inline void Foam::token::operator=(const string& s)
{
    operator=(new string(s));
}

inline void Foam::token::operator=(const label l)
{
    clear();
    type_ = LABEL;
    labelToken_ = l;
}

inline void Foam::token::operator=(const floatScalar s)
{
    clear();
    type_ = FLOAT_SCALAR;
    floatScalarToken_ = s;
}

inline void Foam::token::operator=(const doubleScalar s)
{
    clear();
    type_ = DOUBLE_SCALAR;
    doubleScalarToken_ = s;
}

inline void Foam::token::operator=(const longDoubleScalar s)
{
    clear();
    type_ = LONG_DOUBLE_SCALAR;
    longDoubleScalarTokenPtr_ = new longDoubleScalar(s);
}

inline void Foam::token::operator=(Foam::token::compound* cPtr)
{
    clear();
    type_ = COMPOUND;
    compoundTokenPtr_ = cPtr;
}


inline bool Foam::token::operator==(const token& t) const
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

        case LONG_DOUBLE_SCALAR:
            return equal
            (
                *longDoubleScalarTokenPtr_,
                *t.longDoubleScalarTokenPtr_
            );

        case COMPOUND:
            return compoundTokenPtr_ == t.compoundTokenPtr_;

        case token::ERROR:
            return true;
    }

    return false;
}

inline bool Foam::token::operator==(const punctuationToken p) const
{
    return (type_ == PUNCTUATION && punctuationToken_ == p);
}

inline bool Foam::token::operator==(const word& w) const
{
    return (type_ == WORD && wordToken() == w);
}

inline bool Foam::token::operator==(const string& s) const
{
    return
    (
        (type_ == STRING || type_ == VARIABLE || type_ == VERBATIMSTRING)
     && stringToken() == s
    );
}

inline bool Foam::token::operator==(const label l) const
{
    return (type_ == LABEL && labelToken_ == l);
}

inline bool Foam::token::operator==(const floatScalar s) const
{
    return (type_ == FLOAT_SCALAR && equal(floatScalarToken_, s));
}

inline bool Foam::token::operator==(const doubleScalar s) const
{
    return (type_ == DOUBLE_SCALAR && equal(doubleScalarToken_, s));
}

inline bool Foam::token::operator==(const longDoubleScalar s) const
{
    return
    (
        type_ == LONG_DOUBLE_SCALAR && equal(*longDoubleScalarTokenPtr_, s)
    );
}

inline bool Foam::token::operator!=(const token& t) const
{
    return !operator==(t);
}

inline bool Foam::token::operator!=(const punctuationToken p) const
{
    return !operator==(p);
}

inline bool Foam::token::operator!=(const word& w) const
{
    return !operator==(w);
}

inline bool Foam::token::operator!=(const string& s) const
{
    return !operator==(s);
}

inline bool Foam::token::operator!=(const floatScalar s) const
{
    return !operator==(s);
}

inline bool Foam::token::operator!=(const doubleScalar s) const
{
    return !operator==(s);
}

inline bool Foam::token::operator!=(const longDoubleScalar s) const
{
    return !operator==(s);
}

inline bool Foam::token::operator!=(const label l) const
{
    return !operator==(l);
}


// ************************************************************************* //
