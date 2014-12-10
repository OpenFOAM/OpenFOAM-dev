/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dimensionSet, 1);
const scalar dimensionSet::smallExponent = SMALL;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dimensionSet::dimensionSet
(
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles,
    const scalar current,
    const scalar luminousIntensity
)
{
    exponents_[MASS] = mass;
    exponents_[LENGTH] = length;
    exponents_[TIME] = time;
    exponents_[TEMPERATURE] = temperature;
    exponents_[MOLES] = moles;
    exponents_[CURRENT] = current;
    exponents_[LUMINOUS_INTENSITY] = luminousIntensity;
}


Foam::dimensionSet::dimensionSet
(
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles
)
{
    exponents_[MASS] = mass;
    exponents_[LENGTH] = length;
    exponents_[TIME] = time;
    exponents_[TEMPERATURE] = temperature;
    exponents_[MOLES] = moles;
    exponents_[CURRENT] = 0;
    exponents_[LUMINOUS_INTENSITY] = 0;
}


Foam::dimensionSet::dimensionSet(const dimensionSet& ds)
{
    reset(ds);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dimensionSet::dimensionless() const
{
    for (int Dimension=0; Dimension<nDimensions; ++Dimension)
    {
        // ie, mag(exponents_[Dimension]) > smallExponent
        if
        (
            exponents_[Dimension] > smallExponent
         || exponents_[Dimension] < -smallExponent
        )
        {
            return false;
        }
    }

    return true;
}


void Foam::dimensionSet::reset(const dimensionSet& ds)
{
    for (int Dimension=0; Dimension<nDimensions; ++Dimension)
    {
        exponents_[Dimension] = ds.exponents_[Dimension];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::scalar Foam::dimensionSet::operator[](const dimensionType type) const
{
    return exponents_[type];
}


Foam::scalar& Foam::dimensionSet::operator[](const dimensionType type)
{
    return exponents_[type];
}


Foam::scalar Foam::dimensionSet::operator[](const label type) const
{
    return exponents_[type];
}


Foam::scalar& Foam::dimensionSet::operator[](const label type)
{
    return exponents_[type];
}


bool Foam::dimensionSet::operator==(const dimensionSet& ds) const
{
    for (int Dimension=0; Dimension < nDimensions; ++Dimension)
    {
        if
        (
            mag(exponents_[Dimension] - ds.exponents_[Dimension])
          > smallExponent
        )
        {
            return false;
        }
    }

    return true;
}


bool Foam::dimensionSet::operator!=(const dimensionSet& ds) const
{
    return !(operator==(ds));
}


bool Foam::dimensionSet::operator=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorIn("dimensionSet::operator=(const dimensionSet&) const")
            << "Different dimensions for =" << endl
            << "     dimensions : " << *this << " = " << ds << endl
            << abort(FatalError);
    }

    return true;
}


bool Foam::dimensionSet::operator+=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorIn("dimensionSet::operator+=(const dimensionSet&) const")
            << "Different dimensions for +=" << endl
            << "     dimensions : " << *this << " = " << ds << endl
            << abort(FatalError);
    }

    return true;
}


bool Foam::dimensionSet::operator-=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorIn("dimensionSet::operator-=(const dimensionSet&) const")
            << "Different dimensions for -=" << endl
            << "     dimensions : " << *this << " = " << ds << endl
            << abort(FatalError);
    }

    return true;
}


bool Foam::dimensionSet::operator*=(const dimensionSet& ds)
{
    reset((*this)*ds);

    return true;
}


bool Foam::dimensionSet::operator/=(const dimensionSet& ds)
{
    reset((*this)/ds);

    return true;
}


// * * * * * * * * * * * * * * * Friend functions * * * * * * * * * * * * * * //

Foam::dimensionSet Foam::max(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn("max(const dimensionSet&, const dimensionSet&)")
            << "Arguments of max have different dimensions" << endl
            << "     dimensions : " << ds1 << " and " << ds2 << endl
            << abort(FatalError);
    }

    return ds1;
}


Foam::dimensionSet Foam::min(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn("min(const dimensionSet&, const dimensionSet&)")
            << "Arguments of min have different dimensions" << endl
            << "     dimensions : " << ds1 << " and " << ds2 << endl
            << abort(FatalError);
    }

    return ds1;
}


Foam::dimensionSet Foam::cmptMultiply
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


Foam::dimensionSet Foam::cmptDivide
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1/ds2;
}


Foam::dimensionSet Foam::pow(const dimensionSet& ds, const scalar p)
{
    dimensionSet dimPow
    (
        ds[dimensionSet::MASS]*p,
        ds[dimensionSet::LENGTH]*p,
        ds[dimensionSet::TIME]*p,
        ds[dimensionSet::TEMPERATURE]*p,
        ds[dimensionSet::MOLES]*p,
        ds[dimensionSet::CURRENT]*p,
        ds[dimensionSet::LUMINOUS_INTENSITY]*p
    );

    return dimPow;
}


Foam::dimensionSet Foam::pow
(
    const dimensionSet& ds,
    const dimensionedScalar& dS
)
{
    if (dimensionSet::debug && !dS.dimensions().dimensionless())
    {
        FatalErrorIn("pow(const dimensionSet&, const dimensionedScalar&)")
            << "Exponent of pow is not dimensionless"
            << abort(FatalError);
    }

    dimensionSet dimPow
    (
        ds[dimensionSet::MASS]*dS.value(),
        ds[dimensionSet::LENGTH]*dS.value(),
        ds[dimensionSet::TIME]*dS.value(),
        ds[dimensionSet::TEMPERATURE]*dS.value(),
        ds[dimensionSet::MOLES]*dS.value(),
        ds[dimensionSet::CURRENT]*dS.value(),
        ds[dimensionSet::LUMINOUS_INTENSITY]*dS.value()
    );

    return dimPow;
}


Foam::dimensionSet Foam::pow
(
    const dimensionedScalar& dS,
    const dimensionSet& ds
)
{
    if
    (
        dimensionSet::debug
     && !dS.dimensions().dimensionless()
     && !ds.dimensionless()
    )
    {
        FatalErrorIn("pow(const dimensionedScalar&, const dimensionSet&)")
            << "Argument or exponent of pow not dimensionless" << endl
            << abort(FatalError);
    }

    return ds;
}


Foam::dimensionSet Foam::sqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}


Foam::dimensionSet Foam::pow3(const dimensionSet& ds)
{
    return pow(ds, 3);
}


Foam::dimensionSet Foam::pow4(const dimensionSet& ds)
{
    return pow(ds, 4);
}


Foam::dimensionSet Foam::pow5(const dimensionSet& ds)
{
    return pow(ds, 5);
}


Foam::dimensionSet Foam::pow6(const dimensionSet& ds)
{
    return pow(ds, 6);
}


Foam::dimensionSet Foam::pow025(const dimensionSet& ds)
{
    return sqrt(sqrt(ds));
}


Foam::dimensionSet Foam::sqrt(const dimensionSet& ds)
{
    return pow(ds, 0.5);
}


Foam::dimensionSet Foam::cbrt(const dimensionSet& ds)
{
    return pow(ds, 1.0/3.0);
}


Foam::dimensionSet Foam::magSqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}


Foam::dimensionSet Foam::mag(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::sign(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::pos(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::neg(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::inv(const dimensionSet& ds)
{
    return dimless/ds;
}


Foam::dimensionSet Foam::trans(const dimensionSet& ds)
{
    if (dimensionSet::debug && !ds.dimensionless())
    {
        FatalErrorIn("trans(const dimensionSet&)")
            << "Argument of trancendental function not dimensionless"
            << abort(FatalError);
    }

    return ds;
}


Foam::dimensionSet Foam::atan2(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn("atan2(const dimensionSet&, const dimensionSet&)")
            << "Arguments of atan2 have different dimensions" << endl
            << "     dimensions : " << ds1 << " and " << ds2 << endl
            << abort(FatalError);
    }

    return dimless;
}


Foam::dimensionSet Foam::transform(const dimensionSet& ds)
{
    return ds;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::dimensionSet Foam::operator-(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::operator+
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet dimSum(ds1);

    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn
            ("operator+(const dimensionSet&, const dimensionSet&)")
            << "LHS and RHS of + have different dimensions" << endl
            << "     dimensions : " << ds1 << " + " << ds2 << endl
            << abort(FatalError);
    }

    return dimSum;
}


Foam::dimensionSet Foam::operator-
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet dimDifference(ds1);

    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn
            ("operator-(const dimensionSet&, const dimensionSet&)")
            << "LHS and RHS of - have different dimensions" << endl
            << "     dimensions : " << ds1 << " - " << ds2 << endl
            << abort(FatalError);
    }

    return dimDifference;
}


Foam::dimensionSet Foam::operator*
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet dimProduct(ds1);

    for (int Dimension=0; Dimension<dimensionSet::nDimensions; Dimension++)
    {
        dimProduct.exponents_[Dimension] += ds2.exponents_[Dimension];
    }

    return dimProduct;
}


Foam::dimensionSet Foam::operator/
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet dimQuotient(ds1);

    for (int Dimension=0; Dimension<dimensionSet::nDimensions; Dimension++)
    {
        dimQuotient.exponents_[Dimension] -= ds2.exponents_[Dimension];
    }

    return dimQuotient;
}


Foam::dimensionSet Foam::operator&
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


Foam::dimensionSet Foam::operator^
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


Foam::dimensionSet Foam::operator&&
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


// ************************************************************************* //
