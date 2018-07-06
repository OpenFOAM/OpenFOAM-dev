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

#include "polynomialFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polynomialFunction, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


Foam::polynomialFunction Foam::polynomialFunction::cloneIntegral
(
    const polynomialFunction& poly,
    const scalar intConstant
)
{
    polynomialFunction newPoly(poly.size()+1);

    newPoly[0] = intConstant;
    forAll(poly, i)
    {
        newPoly[i+1] = poly[i]/(i + 1);
    }

    return newPoly;
}


Foam::polynomialFunction Foam::polynomialFunction::cloneIntegralMinus1
(
    const polynomialFunction& poly,
    const scalar intConstant
)
{
    polynomialFunction newPoly(poly.size()+1);

    if (poly[0] > vSmall)
    {
        newPoly.logActive_ = true;
        newPoly.logCoeff_  = poly[0];
    }

    newPoly[0] = intConstant;
    for (label i=1; i < poly.size(); ++i)
    {
        newPoly[i] = poly[i]/i;
    }

    return newPoly;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polynomialFunction::polynomialFunction(const label order)
:
    scalarList(order, 0.0),
    logActive_(false),
    logCoeff_(0.0)
{
    if (this->empty())
    {
        FatalErrorInFunction
            << "polynomialFunction coefficients are invalid (empty)"
            << nl << exit(FatalError);
    }
}


Foam::polynomialFunction::polynomialFunction(const polynomialFunction& poly)
:
    scalarList(poly),
    logActive_(poly.logActive_),
    logCoeff_(poly.logCoeff_)
{}


Foam::polynomialFunction::polynomialFunction(const UList<scalar>& coeffs)
:
    scalarList(coeffs),
    logActive_(false),
    logCoeff_(0.0)
{
    if (this->empty())
    {
        FatalErrorInFunction
            << "polynomialFunction coefficients are invalid (empty)"
            << nl << exit(FatalError);
    }
}


Foam::polynomialFunction::polynomialFunction(Istream& is)
:
    scalarList(is),
    logActive_(false),
    logCoeff_(0.0)
{
    if (this->empty())
    {
        FatalErrorInFunction
            << "polynomialFunction coefficients are invalid (empty)"
            << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polynomialFunction::~polynomialFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polynomialFunction::logActive() const
{
    return logActive_;
}


Foam::scalar Foam::polynomialFunction::logCoeff() const
{
    return logCoeff_;
}


Foam::scalar Foam::polynomialFunction::value(const scalar x) const
{
    const scalarList& coeffs = *this;
    scalar val = coeffs[0];

    // avoid costly pow() in calculation
    scalar powX = x;
    for (label i=1; i<coeffs.size(); ++i)
    {
        val += coeffs[i]*powX;
        powX *= x;
    }

    if (logActive_)
    {
        val += this->logCoeff_*log(x);
    }

    return val;
}


Foam::scalar Foam::polynomialFunction::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    const scalarList& coeffs = *this;

    if (logActive_)
    {
        FatalErrorInFunction
            << "Cannot integrate polynomial with logarithmic coefficients"
            << nl << abort(FatalError);
    }

    // avoid costly pow() in calculation
    scalar powX1 = x1;
    scalar powX2 = x2;

    scalar val = coeffs[0]*(powX2 - powX1);
    for (label i=1; i<coeffs.size(); ++i)
    {
        val += coeffs[i]/(i + 1)*(powX2 - powX1);
        powX1 *= x1;
        powX2 *= x2;
    }

    return val;
}


Foam::polynomialFunction
Foam::polynomialFunction::integral(const scalar intConstant) const
{
    return cloneIntegral(*this, intConstant);
}


Foam::polynomialFunction
Foam::polynomialFunction::integralMinus1(const scalar intConstant) const
{
    return cloneIntegralMinus1(*this, intConstant);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::polynomialFunction&
Foam::polynomialFunction::operator+=(const polynomialFunction& poly)
{
    scalarList& coeffs = *this;

    if (coeffs.size() > poly.size())
    {
        forAll(poly, i)
        {
            coeffs[i] += poly[i];
        }
    }
    else
    {
        coeffs.setSize(poly.size(), 0.0);

        forAll(coeffs, i)
        {
            coeffs[i] += poly[i];
        }
    }

    return *this;
}


Foam::polynomialFunction&
Foam::polynomialFunction::operator-=(const polynomialFunction& poly)
{
    scalarList& coeffs = *this;

    if (coeffs.size() > poly.size())
    {
        forAll(poly, i)
        {
            coeffs[i] -= poly[i];
        }
    }
    else
    {
        coeffs.setSize(poly.size(), 0.0);

        forAll(coeffs, i)
        {
            coeffs[i] -= poly[i];
        }
    }

    return *this;
}


Foam::polynomialFunction&
Foam::polynomialFunction::operator*=(const scalar s)
{
    scalarList& coeffs = *this;
    forAll(coeffs, i)
    {
        coeffs[i] *= s;
    }

    return *this;
}


Foam::polynomialFunction&
Foam::polynomialFunction::operator/=(const scalar s)
{
    scalarList& coeffs = *this;
    forAll(coeffs, i)
    {
        coeffs[i] /= s;
    }

    return *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polynomialFunction& poly)
{
    // output like VectorSpace
    os  << token::BEGIN_LIST;

    if (!poly.empty())
    {
        for (int i=0; i<poly.size()-1; i++)
        {
            os  << poly[i] << token::SPACE;
        }
        os  << poly.last();
    }
    os  << token::END_LIST;


    // Check state of Ostream
    os.check("operator<<(Ostream&, const polynomialFunction&)");

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::polynomialFunction
Foam::operator+
(
    const polynomialFunction& p1,
    const polynomialFunction& p2
)
{
    polynomialFunction poly(p1);
    return poly += p2;
}


Foam::polynomialFunction
Foam::operator-
(
    const polynomialFunction& p1,
    const polynomialFunction& p2
)
{
    polynomialFunction poly(p1);
    return poly -= p2;
}


Foam::polynomialFunction
Foam::operator*
(
    const scalar s,
    const polynomialFunction& p
)
{
    polynomialFunction poly(p);
    return poly *= s;
}


Foam::polynomialFunction
Foam::operator/
(
    const scalar s,
    const polynomialFunction& p
)
{
    polynomialFunction poly(p);
    return poly /= s;
}


Foam::polynomialFunction
Foam::operator*
(
    const polynomialFunction& p,
    const scalar s
)
{
    polynomialFunction poly(p);
    return poly *= s;
}


Foam::polynomialFunction
Foam::operator/
(
    const polynomialFunction& p,
    const scalar s
)
{
    polynomialFunction poly(p);
    return poly /= s;
}


// ************************************************************************* //
