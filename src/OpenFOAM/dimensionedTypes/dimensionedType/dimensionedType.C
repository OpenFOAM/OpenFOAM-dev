/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "dimensionedType.H"
#include "pTraits.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::dimensioned<Type>::initialise(Istream& is)
{
    token nextToken(is);
    is.putBack(nextToken);

    // Check if the original format is used in which the name is provided
    // and reset the name to that read
    if (nextToken.isWord())
    {
        is >> name_;
        is >> nextToken;
        is.putBack(nextToken);
    }

    // If the dimensions are provided compare with the argument
    scalar multiplier = 1.0;

    if (nextToken == token::BEGIN_SQR)
    {
        dimensionSet dims(dimless);
        dims.read(is, multiplier);

        if (dims != dimensions_)
        {
            FatalIOErrorInFunction(is)
                << "The dimensions " << dims
                << " provided do not match the required dimensions "
                << dimensions_ << abort(FatalIOError);
        }
    }

    is >> value_;
    value_ *= multiplier;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dimSet,
    const Type& t
)
:
    name_(name),
    dimensions_(dimSet),
    value_(t)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const dimensionSet& dimSet,
    const Type& t
)
:
    name_(::Foam::name(t)),
    dimensions_(dimSet),
    value_(t)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned(const Type& t)
:
    name_(::Foam::name(t)),
    dimensions_(dimless),
    value_(t)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensioned<Type>& dt
)
:
    name_(name),
    dimensions_(dt.dimensions_),
    value_(dt.value_)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    Istream& is
)
:
    dimensions_(dimless)
{
    read(is);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    Istream& is
)
:
    name_(name),
    dimensions_(dimless)
{
    scalar multiplier;
    dimensions_.read(is, multiplier);
    is >> value_;
    value_ *= multiplier;
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dimSet,
    Istream& is
)
:
    name_(name),
    dimensions_(dimSet),
    value_(Zero)
{
    initialise(is);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dimSet,
    const dictionary& dict
)
:
    name_(name),
    dimensions_(dimSet),
    value_(Zero)
{
    initialise(dict.lookup(name));
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
()
:
    name_("undefined"),
    dimensions_(dimless),
    value_(Zero)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::lookupOrDefault
(
    const word& name,
    const dictionary& dict,
    const dimensionSet& dims,
    const Type& defaultValue
)
{
    if (dict.found(name))
    {
        return dimensioned<Type>(name, dims, dict.lookup(name));
    }
    else
    {
        return dimensioned<Type>(name, dims, defaultValue);
    }
}


template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::lookupOrDefault
(
    const word& name,
    const dictionary& dict,
    const Type& defaultValue
)
{
    return lookupOrDefault(name, dict, dimless, defaultValue);
}


template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const dimensionSet& dims,
    const Type& defaultValue
)
{
    Type value = dict.lookupOrAddDefault<Type>(name, defaultValue);
    return dimensioned<Type>(name, dims, value);
}


template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Type& defaultValue
)
{
    return lookupOrAddToDict(name, dict, dimless, defaultValue);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::dimensioned<Type>::name() const
{
    return name_;
}

template<class Type>
Foam::word& Foam::dimensioned<Type>::name()
{
    return name_;
}


template<class Type>
const Foam::dimensionSet& Foam::dimensioned<Type>::dimensions() const
{
    return dimensions_;
}

template<class Type>
Foam::dimensionSet& Foam::dimensioned<Type>::dimensions()
{
    return dimensions_;
}


template<class Type>
const Type& Foam::dimensioned<Type>::value() const
{
    return value_;
}

template<class Type>
Type& Foam::dimensioned<Type>::value()
{
    return value_;
}


template<class Type>
Foam::dimensioned<typename Foam::dimensioned<Type>::cmptType>
Foam::dimensioned<Type>::component
(
    const direction d
) const
{
    return dimensioned<cmptType>
    (
        name_ + ".component(" + Foam::name(d) + ')',
        dimensions_,
        value_.component(d)
    );
}


template<class Type>
void Foam::dimensioned<Type>::replace
(
    const direction d,
    const dimensioned<typename dimensioned<Type>::cmptType>& dc
)
{
    dimensions_ = dc.dimensions();
    value_.replace(d, dc.value());
}


template<class Type>
void Foam::dimensioned<Type>::read(const dictionary& dict)
{
    initialise(dict.lookup(name_));
}


template<class Type>
bool Foam::dimensioned<Type>::readIfPresent(const dictionary& dict)
{
    return dict.readIfPresent(name_, value_);
}


template<class Type>
Foam::Istream& Foam::dimensioned<Type>::read(Istream& is)
{
    // If the name is present, read it
    token nextToken(is);
    is.putBack(nextToken);
    if (nextToken.isWord())
    {
        is >> name_;
    }

    // Read the dimensions and multiplier
    scalar multiplier;
    dimensions_.read(is, multiplier);

    // Read and scale the value
    is >> value_;
    value_ *= multiplier;

    // If the name is not present, set it
    if (!nextToken.isWord())
    {
        name_ = Foam::name(value_);
    }

    // Check state of Istream
    is.check
    (
        "Istream& dimensioned<Type>::read(Istream& is)"
    );

    return is;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<typename Foam::dimensioned<Type>::cmptType>
Foam::dimensioned<Type>::operator[]
(
    const direction d
) const
{
    return component(d);
}


template<class Type>
void Foam::dimensioned<Type>::operator+=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ += dt.dimensions_;
    value_ += dt.value_;
}


template<class Type>
void Foam::dimensioned<Type>::operator-=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ -= dt.dimensions_;
    value_ -= dt.value_;
}


template<class Type>
void Foam::dimensioned<Type>::operator*=
(
    const scalar s
)
{
    value_ *= s;
}


template<class Type>
void Foam::dimensioned<Type>::operator/=
(
    const scalar s
)
{
    value_ /= s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type, Foam::direction r>
Foam::dimensioned<typename Foam::powProduct<Type, r>::type>
Foam::pow(const dimensioned<Type>& dt, typename powProduct<Type, r>::type)
{
    return dimensioned<typename powProduct<Type, r>::type>
    (
        "pow(" + dt.name() + ',' + name(r) + ')',
        pow(dt.dimensions(), r),
        pow(dt.value(), 2)
    );
}


template<class Type>
Foam::dimensioned<typename Foam::outerProduct<Type, Type>::type>
Foam::sqr(const dimensioned<Type>& dt)
{
    return dimensioned<typename outerProduct<Type, Type>::type>
    (
        "sqr(" + dt.name() + ')',
        sqr(dt.dimensions()),
        sqr(dt.value())
    );
}

template<class Type>
Foam::dimensioned<Foam::scalar> Foam::magSqr(const dimensioned<Type>& dt)
{
    return dimensioned<scalar>
    (
        "magSqr(" + dt.name() + ')',
        magSqr(dt.dimensions()),
        magSqr(dt.value())
    );
}

template<class Type>
Foam::dimensioned<Foam::scalar> Foam::mag(const dimensioned<Type>& dt)
{
    return dimensioned<scalar>
    (
        "mag(" + dt.name() + ')',
        dt.dimensions(),
        mag(dt.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::cmptMultiply
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        "cmptMultiply(" + dt1.name() + ',' + dt2.name() + ')',
        cmptMultiply(dt1.dimensions(), dt2.dimensions()),
        cmptMultiply(dt1.value(), dt2.value())
    );
}

template<class Type>
Foam::dimensioned<Type> Foam::cmptDivide
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        "cmptDivide(" + dt1.name() + ',' + dt2.name() + ')',
        cmptDivide(dt1.dimensions(), dt2.dimensions()),
        cmptDivide(dt1.value(), dt2.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::max
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    if (dt1.dimensions() != dt2.dimensions())
    {
        FatalErrorInFunction
            << "dimensions of arguments are not equal"
            << abort(FatalError);
    }

    return dimensioned<Type>
    (
        "max(" + dt1.name() + ',' + dt2.name() + ')',
        dt1.dimensions(),
        max(dt1.value(), dt2.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::min
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    if (dt1.dimensions() != dt2.dimensions())
    {
        FatalErrorInFunction
            << "dimensions of arguments are not equal"
            << abort(FatalError);
    }

    return dimensioned<Type>
    (
        "min(" + dt1.name() + ',' + dt2.name() + ')',
        dt1.dimensions(),
        min(dt1.value(), dt2.value())
    );
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::writeEntry(Ostream& os, const dimensioned<Type>& dt)
{
    os << dt;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, dimensioned<Type>& dt)
{
    token nextToken(is);
    is.putBack(nextToken);

    // Check if the original format is used in which the name is provided
    // and reset the name to that read
    if (nextToken.isWord())
    {
        is >> dt.name_;
        is >> nextToken;
        is.putBack(nextToken);
    }

    // If the dimensions are provided reset the dimensions to those read
    scalar multiplier = 1.0;
    if (nextToken == token::BEGIN_SQR)
    {
        dt.dimensions_.read(is, multiplier);
    }

    // Read the value
    is >> dt.value_;
    dt.value_ *= multiplier;

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, dimensioned<Type>&)");

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const dimensioned<Type>& dt)
{
    // Write the name
    os << dt.name() << token::SPACE;

    // Write the dimensions
    dt.dimensions().write(os);

    os << token::SPACE;

    // Write the value
    os << dt.value();

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const dimensioned<Type>&)");

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
bool Foam::operator>
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dt1.value() > dt2.value();
}


template<class Type>
bool Foam::operator<
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dt1.value() < dt2.value();
}


template<class Type>
Foam::dimensioned<Type> Foam::operator+
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        '(' + dt1.name() + '+' + dt2.name() + ')',
        dt1.dimensions() + dt2.dimensions(),
        dt1.value() + dt2.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator-(const dimensioned<Type>& dt)
{
    return dimensioned<Type>
    (
        '-' + dt.name(),
        dt.dimensions(),
        -dt.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator-
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        '(' + dt1.name() + '-' + dt2.name() + ')',
        dt1.dimensions() - dt2.dimensions(),
        dt1.value() - dt2.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const dimensioned<Type>& dt
)
{
    return dimensioned<Type>
    (
        '(' + ds.name() + '*' + dt.name() + ')',
        ds.dimensions() * dt.dimensions(),
        ds.value() * dt.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator/
(
    const dimensioned<Type>& dt,
    const dimensioned<scalar>& ds
)
{
    return dimensioned<Type>
    (
        '(' + dt.name() + '|' + ds.name() + ')',
        dt.dimensions()/ds.dimensions(),
        dt.value()/ds.value()
    );
}


#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
Foam::dimensioned<typename Foam::product<Type1, Type2>::type>                  \
Foam::operator op                                                              \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    return dimensioned<typename product<Type1, Type2>::type>                   \
    (                                                                          \
        '(' + dt1.name() + #op + dt2.name() + ')',                             \
        dt1.dimensions() op dt2.dimensions(),                                  \
        dt1.value() op dt2.value()                                             \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, Foam::direction nCmpt>            \
Foam::dimensioned<typename Foam::product<Type, Form>::type>                    \
Foam::operator op                                                              \
(                                                                              \
    const dimensioned<Type>& dt1,                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& t2                                     \
)                                                                              \
{                                                                              \
    return dimensioned<typename product<Type, Form>::type>                     \
    (                                                                          \
        '(' + dt1.name() + #op + name(t2) + ')',                               \
        dt1.dimensions(),                                                      \
        dt1.value() op static_cast<const Form&>(t2)                            \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, Foam::direction nCmpt>            \
Foam::dimensioned<typename Foam::product<Form, Type>::type>                    \
Foam::operator op                                                              \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& t1,                                    \
    const dimensioned<Type>& dt2                                               \
)                                                                              \
{                                                                              \
    return dimensioned<typename product<Form, Type>::type>                     \
    (                                                                          \
        '(' + name(t1) + #op + dt2.name() + ')',                               \
        dt2.dimensions(),                                                      \
        static_cast<const Form&>(t1) op dt2.value()                            \
    );                                                                         \
}


PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// ************************************************************************* //
