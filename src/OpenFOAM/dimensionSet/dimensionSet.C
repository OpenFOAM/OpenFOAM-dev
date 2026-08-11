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

#include "dimensionSet.H"
#include "dimensions.H"
#include "dimensionedScalar.H"
#include "NamedEnum.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dimensionSet, 1);
}

const Foam::autoPtr<Foam::NamedEnum<Foam::dimensionSet::dimensionType, 7>>
    dimensionTypeNamesPtr_
    (
        new Foam::NamedEnum<Foam::dimensionSet::dimensionType, 7>
        {
            "mass",
            "length",
            "time",
            "temperature",
            "moles",
            "current",
            "luminousIntensity"
        }
    );

const Foam::NamedEnum<Foam::dimensionSet::dimensionType, 7>&
    Foam::dimensionSet::dimensionTypeNames_ = dimensionTypeNamesPtr_();

const Foam::scalar Foam::dimensionSet::smallExponent = rootSmall;


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
    const word& name,
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles,
    const scalar current,
    const scalar luminousIntensity
)
:
    dimensionSet
    (
        mass,
        length,
        time,
        temperature,
        moles,
        current,
        luminousIntensity
    )
{
    dimensions::table.insert(name, *this);
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


Foam::dimensionSet::dimensionSet(const word& name, const dimensionSet& ds)
{
    reset(ds);
    dimensions::table.insert(name, *this);
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


// * * * * * * * * * * * * * * * Friend functions * * * * * * * * * * * * * * //

Foam::dimensionSet Foam::pow
(
    const dimensionSet& ds,
    const dimensionedScalar& dS
)
{
    if (dimensionSet::debug && !dS.dimensions().dimensionless())
    {
        FatalErrorInFunction
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
        FatalErrorInFunction
            << "Argument or exponent of pow not dimensionless" << endl
            << abort(FatalError);
    }

    return ds;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

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


// ************************************************************************* //
