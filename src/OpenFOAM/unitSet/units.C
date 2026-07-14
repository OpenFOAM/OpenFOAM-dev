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

#include "units.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "symbols.H"
#include "mathematicalConstants.H"

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
                debug::configDict().found("units")
              ? "units"
              : debug::configDict().found("UnitSets")
              ? "UnitSets"
              : debug::configDict().found("DimensionSets")
              ? "DimensionSets"
              : "units",
                cachedPtr
            )
        );
    }

    return *unitsDictPtr_;
}

HashTable<unitSet>* addedUnitsPtr_(nullptr);
HashTable<unitSet>* unitsPtr_(nullptr);

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
namespace units
{

unitSet newFraction()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 1, 0, 1);
}

unitSet newPercent()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 1, 0, 0.01);
}

unitSet newRadians()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 0, 1, 1);
}

unitSet newRotations()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 0, 1, 2*pi);
}

unitSet newDegrees()
{
    return unitSet(dimensionSet(0, 0, 0, 0, 0), 0, 1, pi/180);
}

unitSet length_(dimensionSet(0, 1, 0, 0, 0), 0, 0, 1);

}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::unitSet Foam::units::unitless(unitSet::newUnitless());

const Foam::unitSet Foam::units::any(unitSet::newAny());
const Foam::unitSet Foam::units::none(unitSet::newNone());

const Foam::unitSet Foam::units::fraction(newFraction());
const Foam::unitSet Foam::units::percent(newPercent());

const Foam::unitSet Foam::units::radians(newRadians());
const Foam::unitSet Foam::units::rotations(newRotations());
const Foam::unitSet Foam::units::degrees(newDegrees());

const Foam::unitSet& Foam::units::length = Foam::units::length_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::HashTable<Foam::unitSet>& Foam::units::table()
{
    if (!unitsPtr_)
    {
        unitsPtr_ = new HashTable<unitSet>();

        unitsPtr_->insert("%", newPercent());

        unitsPtr_->insert("rad", newRadians());
        unitsPtr_->insert("rot", newRotations());
        unitsPtr_->insert("deg", newDegrees());

        // Get the relevant part of the control dictionary
        const dictionary& unitSetDict =
            unitsDict().subDict(unitsDict().lookup<word>("set"));

        // Read fundamental units
        HashTable<label> fundamentalUnits;
        forAllConstIter(dictionary, unitSetDict, iter)
        {
            ITstream& is = iter().stream();

            const dimensionSet units(is);

            // Check that this is unique
            if (fundamentalUnits.found(iter().keyword()))
            {
                FatalIOErrorInFunction(unitsDict())
                    << "Duplicate-fundamental unit specified"
                    << exit(FatalIOError);
            }

            // Get the index of the dimension of this fundamental unit
            label i = -1;
            for (label j = 0; j < dimensionSet::nDimensions; ++ j)
            {
                if (units[j] == 0);
                else if (units[j] == 1 && i == -1) i = j;
                else i = -2;
            }
            if (i < 0)
            {
                FatalIOErrorInFunction(unitsDict())
                    << "Non-fundamental unit specified" << units
                    << exit(FatalIOError);
            }

            // Store the name
            fundamentalUnits.insert(iter().keyword(), i);

            // Add to the table
            unitsPtr_->insert
            (
                iter().keyword(),
                units*unitSet(dimless, 0, 0, scalar(1))
            );
        }

        // Local function to add a unit to the table with checking
        auto addUnit = []
        (
            const word& name,
            const unitSet& units,
            const scalar multiplier
        )
        {
            const bool ok =
                unitsPtr_->insert
                (
                    name,
                    units*unitSet(dimless, 0, 0, multiplier)
                );

            if (!ok)
            {
                FatalIOErrorInFunction(unitsDict())
                    << "Duplicate unit " << name
                    << " read from dictionary"
                    << exit(FatalIOError);
            }
        };

        // Read conversions to fundamental units. These are defined the wrong
        // way around relative to the chosen unit set, so we have to parse them
        // and then add the reverse conversion to the table.
        forAllConstIter(HashTable<label>, fundamentalUnits, iter)
        {
            if (!unitsDict().found(iter.key())) continue;

            auto readPunctuation = [](symbols::tokeniser& tis, const char p)
            {
                token t(tis.nextToken());
                if (!t.isPunctuation() || t.pToken() != p)
                {
                    FatalIOErrorInFunction(tis.stream())
                        << "Illegal token " << t
                        << ". Expected character '" << p << "'."
                        << exit(FatalIOError);
                }
            };

            ITstream& is = unitsDict().lookup(iter.key());
            symbols::tokeniser tis(is);

            // Read the name and the multiplier, in whatever order
            token t(tis.nextToken());
            tis.putBack(t);
            word derivedUnitName;
            scalar multiplier;
            if (!t.isNumber())
            {
                readPunctuation(tis, token::BEGIN_SQR);
                derivedUnitName = tis.nextToken().wordToken();
                readPunctuation(tis, token::END_SQR);
                multiplier = tis.nextToken().number();
            }
            else
            {
                multiplier = tis.nextToken().number();
                readPunctuation(tis, token::BEGIN_SQR);
                derivedUnitName = tis.nextToken().wordToken();
                readPunctuation(tis, token::END_SQR);
            }

            // Now we have the name and the multiplier, we can add the
            // inverse conversion to the table; i.e., flip the names, and
            // use the reciprocal of the multiplier
            addUnit
            (
                derivedUnitName,
                unitsPtr_->operator[](iter.key()),
                1/multiplier
            );
        }

        // Read other units
        forAllConstIter(dictionary, unitsDict(), iter)
        {
            if
            (
                iter().isDict()
             || iter().keyword() == "set"
             || iter().keyword() == "unitSet"
             || fundamentalUnits.found(iter().keyword())
            ) continue;

            // Read the units and the multiplier, in whatever order
            ITstream& is = iter().stream();
            token t(is);
            is.putBack(t);
            unitSet units(units::any);
            scalar multiplier;
            if (!t.isNumber())
            {
                is >> units;
                multiplier = pTraits<scalar>(is);
            }
            else
            {
                multiplier = pTraits<scalar>(is);
                is >> units;
            }

            // Add the conversion to the table
            addUnit(iter().keyword(), units, multiplier);
        }

        // Add programmatically defined units
        if (addedUnitsPtr_)
        {
            forAllConstIter(HashTable<unitSet>, *addedUnitsPtr_, iter)
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


void Foam::units::add(const word& name, const unitSet& units)
{
    deleteDemandDrivenData(unitsDictPtr_);

    if (!addedUnitsPtr_)
    {
        addedUnitsPtr_ = new HashTable<unitSet>();
    }

    addedUnitsPtr_->insert(name, units);

    deleteDemandDrivenData(unitsPtr_);
}


const Foam::unitSet& Foam::units::lookup(const word& unitName)
{
    return table()[unitName];
}


void Foam::units::setLength(const unitSet& length)
{
    length_.reset(length);
}


void Foam::units::setLength(const scalar length)
{
    length_.reset(unitSet(dimensions::length, 0, 0, length));
}


Foam::scalar Foam::degToRad(const scalar deg)
{
    return units::degrees.toStandard(deg);
}


Foam::scalar Foam::radToDeg(const scalar rad)
{
    return units::degrees.toUser(rad);
}


// ************************************************************************* //
