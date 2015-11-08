/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "PolynomialEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PolynomialEntry<Type>::PolynomialEntry
(
    const word& entryName,
    const dictionary& dict
)
:
    DataEntry<Type>(entryName),
    coeffs_(),
    canIntegrate_(true),
    dimensions_(dimless)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);

    token firstToken(is);
    is.putBack(firstToken);
    if (firstToken == token::BEGIN_SQR)
    {
        is  >> this->dimensions_;
    }

    is  >> coeffs_;

    if (!coeffs_.size())
    {
        FatalErrorInFunction
            << "PolynomialEntry coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }

    forAll(coeffs_, i)
    {
        if (mag(coeffs_[i].second() + pTraits<Type>::one) < ROOTVSMALL)
        {
            canIntegrate_ = false;
            break;
        }
    }

    if (debug)
    {
        if (!canIntegrate_)
        {
            WarningInFunction
                << "PolynomialEntry " << this->name_ << " cannot be integrated"
                << endl;
        }
    }
}


template<class Type>
Foam::PolynomialEntry<Type>::PolynomialEntry
(
    const word& entryName,
    const List<Tuple2<Type, Type> >& coeffs
)
:
    DataEntry<Type>(entryName),
    coeffs_(coeffs),
    canIntegrate_(true),
    dimensions_(dimless)
{
    if (!coeffs_.size())
    {
        FatalErrorInFunction
            << "PolynomialEntry coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }

    forAll(coeffs_, i)
    {
        if (mag(coeffs_[i].second() + 1) < ROOTVSMALL)
        {
            canIntegrate_ = false;
            break;
        }
    }

    if (debug)
    {
        if (!canIntegrate_)
        {
            WarningInFunction
                << "PolynomialEntry " << this->name_ << " cannot be integrated"
                << endl;
        }
    }
}


template<class Type>
Foam::PolynomialEntry<Type>::PolynomialEntry(const PolynomialEntry& poly)
:
    DataEntry<Type>(poly),
    coeffs_(poly.coeffs_),
    canIntegrate_(poly.canIntegrate_),
    dimensions_(poly.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PolynomialEntry<Type>::~PolynomialEntry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PolynomialEntry<Type>::convertTimeBase(const Time& t)
{
    forAll(coeffs_, i)
    {
        Type value = coeffs_[i].first();
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            setComponent(coeffs_[i].first(), cmpt) =
                t.userTimeToTime(component(value, cmpt));
        }
    }
}


template<class Type>
Type Foam::PolynomialEntry<Type>::value(const scalar x) const
{
    Type y(pTraits<Type>::zero);
    forAll(coeffs_, i)
    {
        y += cmptMultiply
        (
            coeffs_[i].first(),
            cmptPow(pTraits<Type>::one*x, coeffs_[i].second())
        );
    }

    return y;
}


template<class Type>
Type Foam::PolynomialEntry<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(pTraits<Type>::zero);

    if (canIntegrate_)
    {
        forAll(coeffs_, i)
        {
            intx += cmptMultiply
            (
                cmptDivide
                (
                    coeffs_[i].first(),
                    coeffs_[i].second() + pTraits<Type>::one
                ),
                cmptPow
                (
                    pTraits<Type>::one*x2,
                    coeffs_[i].second() + pTraits<Type>::one
                )
              - cmptPow
                (
                    pTraits<Type>::one*x1,
                    coeffs_[i].second() + pTraits<Type>::one
                )
            );
        }
    }

    return intx;
}


template<class Type>
Foam::dimensioned<Type> Foam::PolynomialEntry<Type>::dimValue
(
    const scalar x
) const
{
    return dimensioned<Type>("dimensionedValue", dimensions_, value(x));
}


template<class Type>
Foam::dimensioned<Type> Foam::PolynomialEntry<Type>::dimIntegrate
(
    const scalar x1,
    const scalar x2
) const
{
    return dimensioned<Type>
    (
        "dimensionedValue",
        dimensions_,
        integrate(x1, x2)
    );
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "PolynomialEntryIO.C"


// ************************************************************************* //
