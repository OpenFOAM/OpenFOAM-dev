/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "Polynomial1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Polynomial<Type>::Polynomial
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction1<Type, Polynomial<Type>>(name),
    coeffs_(),
    canIntegrate_(true)
{
    if (!dict.found(name))
    {
        dict.lookup("coeffs") >> coeffs_;
    }
    else
    {
        Istream& is(dict.lookup(name));
        word entryType(is);

        if (is.eof())
        {
            dict.lookup("coeffs") >> coeffs_;
        }
        else
        {
            is  >> coeffs_;
        }
    }

    if (!coeffs_.size())
    {
        FatalErrorInFunction
            << "Polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }

    forAll(coeffs_, i)
    {
        if (mag(coeffs_[i].second() + pTraits<Type>::one) < rootVSmall)
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
                << "Polynomial " << this->name_ << " cannot be integrald"
                << endl;
        }
    }
}


template<class Type>
Foam::Function1s::Polynomial<Type>::Polynomial
(
    const word& name,
    const List<Tuple2<Type, Type>>& coeffs
)
:
    FieldFunction1<Type, Polynomial<Type>>(name),
    coeffs_(coeffs),
    canIntegrate_(true)
{
    if (!coeffs_.size())
    {
        FatalErrorInFunction
            << "Polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }

    forAll(coeffs_, i)
    {
        if (mag(coeffs_[i].second() + 1) < rootVSmall)
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
                << "Polynomial " << this->name_ << " cannot be integrald"
                << endl;
        }
    }
}


template<class Type>
Foam::Function1s::Polynomial<Type>::Polynomial(const Polynomial& poly)
:
    FieldFunction1<Type, Polynomial<Type>>(poly),
    coeffs_(poly.coeffs_),
    canIntegrate_(poly.canIntegrate_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Polynomial<Type>::~Polynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::Polynomial<Type>::value(const scalar x) const
{
    Type y(Zero);
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
Type Foam::Function1s::Polynomial<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(Zero);

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
void Foam::Function1s::Polynomial<Type>::write(Ostream& os) const
{
    writeKeyword(os, "coeffs") << coeffs_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
