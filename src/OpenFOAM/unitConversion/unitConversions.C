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

#include "unitConversions.H"
#include "demandDrivenData.H"
#include "dictionary.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

dictionary* unitsDictPtr_(nullptr);

const dictionary& unitsDict()
{
    if (!unitsDictPtr_)
    {
        dictionary* cachedPtr = nullptr;

        unitsDictPtr_ = new dictionary
        (
            debug::switchSet
            (
                debug::controlDict().found("UnitConversions")
              ? "UnitConversions"
              : debug::controlDict().found("DimensionSets")
              ? "DimensionSets"
              : "UnitConversions",
                cachedPtr
            )
        );
    }

    return *unitsDictPtr_;
}

HashTable<unitConversion>* addedUnitsPtr_(nullptr);
HashTable<unitConversion>* unitsPtr_(nullptr);

// Delete the above data at the end of the run
struct deleteUnitsPtr
{
    ~deleteUnitsPtr()
    {
        deleteDemandDrivenData(unitsDictPtr_);
        deleteDemandDrivenData(addedUnitsPtr_);
        deleteDemandDrivenData(unitsPtr_);
    }
};

deleteUnitsPtr deleteUnitsPtr_;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

dimensionSet makeDimless()
{
    return dimensionSet(0, 0, 0, 0, 0);
}

unitConversion makeUnitless()
{
    return unitConversion(makeDimless(), 0, 0, 1);
}

unitConversion makeUnitAny()
{
    return unitConversion(makeDimless(), 0, 0, 0);
}
unitConversion makeUnitNone()
{
    return unitConversion(makeDimless(), 0, 0, -1);
}

unitConversion makeUnitFraction()
{
    return unitConversion(makeDimless(), 1, 0, 1);
}
unitConversion makeUnitPercent()
{
    return unitConversion(makeDimless(), 1, 0, 0.01);
}

unitConversion makeUnitRadians()
{
    return unitConversion(makeDimless(), 0, 1, 1);
}
unitConversion makeUnitRotations()
{
    return unitConversion(makeDimless(), 0, 1, 2*pi);
}
unitConversion makeUnitDegrees()
{
    return unitConversion(makeDimless(), 0, 1, pi/180);
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::unitConversion Foam::unitless(makeUnitless());

const Foam::unitConversion Foam::unitAny(makeUnitAny());
const Foam::unitConversion Foam::unitNone(makeUnitNone());

const Foam::unitConversion Foam::unitFraction(makeUnitFraction());
const Foam::unitConversion Foam::unitPercent(makeUnitPercent());

const Foam::unitConversion Foam::unitRadians(makeUnitRadians());
const Foam::unitConversion Foam::unitRotations(makeUnitRotations());
const Foam::unitConversion Foam::unitDegrees(makeUnitDegrees());


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::degToRad(const scalar deg)
{
    return unitDegrees.toStandard(deg);
}


Foam::scalar Foam::radToDeg(const scalar rad)
{
    return unitDegrees.toUser(rad);
}


void Foam::addUnits(const word& name, const unitConversion& units)
{
    deleteDemandDrivenData(unitsDictPtr_);

    if (!addedUnitsPtr_)
    {
        addedUnitsPtr_ = new HashTable<unitConversion>();
    }

    addedUnitsPtr_->insert(name, units);

    deleteDemandDrivenData(unitsPtr_);
}


const Foam::HashTable<Foam::unitConversion>& Foam::units()
{
    if (!unitsPtr_)
    {
        unitsPtr_ = new HashTable<unitConversion>();

        unitsPtr_->insert("%", makeUnitPercent());

        unitsPtr_->insert("rad", makeUnitRadians());
        unitsPtr_->insert("rot", makeUnitRotations());
        unitsPtr_->insert("deg", makeUnitDegrees());

        // Get the relevant part of the control dictionary
        const dictionary& unitSetDict =
            unitsDict().subDict(unitsDict().lookup<word>("unitSet") + "Coeffs");

        // Add units from the control dictionary
        forAllConstIter(dictionary, unitSetDict, iter)
        {
            ITstream& is = iter().stream();

            const unitConversion units(is);
            const scalar multiplier = pTraits<scalar>(is);

            const bool ok =
                unitsPtr_->insert
                (
                    iter().keyword(),
                    units*unitConversion(dimless, 0, 0, multiplier)
                );

            if (!ok)
            {
                FatalIOErrorInFunction(unitsDict())
                    << "Duplicate unit " << iter().keyword()
                    << " read from dictionary"
                    << exit(FatalIOError);
            }
        }

        // Add programmatically defined units
        if (addedUnitsPtr_)
        {
            forAllConstIter(HashTable<unitConversion>, *addedUnitsPtr_, iter)
            {
                const bool ok = unitsPtr_->insert(iter.key(), iter());

                if (!ok)
                {
                    FatalIOErrorInFunction(unitsDict())
                        << "Duplicate unit " << iter.key()
                        << " added to dictionary"
                        << exit(FatalIOError);
                }
            }
        }
    }

    return *unitsPtr_;
}


// ************************************************************************* //
