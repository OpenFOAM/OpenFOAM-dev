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

const Foam::unitConversion Foam::unitless(dimless, 0, 0, 1);

const Foam::unitConversion Foam::unitAny(dimless, 0, 0, 0);
const Foam::unitConversion Foam::unitNone(dimless, 0, 0, -1);

const Foam::unitConversion Foam::unitFraction(dimless, 1, 0, 1);
const Foam::unitConversion Foam::unitPercent(dimless, 1, 0, 0.01);

const Foam::unitConversion Foam::unitRadians(dimless, 0, 1, 1);
const Foam::unitConversion Foam::unitRotations(dimless, 0, 1, 2*pi);
const Foam::unitConversion Foam::unitDegrees(dimless, 0, 1, pi/180);


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

        unitsPtr_->insert("%", unitPercent);

        unitsPtr_->insert("rad", unitRadians);
        unitsPtr_->insert("rot", unitRotations);
        unitsPtr_->insert("deg", unitDegrees);

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
