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

#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Since dimensionSystems() can be reread we actually store a copy of
// the controlDict subDict (v.s. a reference to the subDict for e.g.
// dimensionedConstants)
dictionary* dimensionSystemsPtr_(nullptr);
HashTable<dimensionedScalar>* unitSetPtr_(nullptr);
dimensionSets* writeUnitSetPtr_(nullptr);

// Delete the above data at the end of the run
struct deleteDimensionSystemsPtr
{
    ~deleteDimensionSystemsPtr()
    {
        deleteDemandDrivenData(dimensionSystemsPtr_);
        deleteDemandDrivenData(unitSetPtr_);
        deleteDemandDrivenData(writeUnitSetPtr_);
    }
};

deleteDimensionSystemsPtr deleteDimensionSystemsPtr_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::dimensionSystems()
{
    if (!dimensionSystemsPtr_)
    {
        dictionary* cachedPtr = nullptr;
        dimensionSystemsPtr_ = new dictionary
        (
            debug::switchSet
            (
                "DimensionSets",
                cachedPtr
            )
        );
    }
    return *dimensionSystemsPtr_;
}


const Foam::HashTable<Foam::dimensionedScalar>& Foam::unitSet()
{
    if (!unitSetPtr_)
    {
        const dictionary& dict = dimensionSystems();

        if (!dict.found("unitSet"))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find unitSet in dictionary " << dict.name()
                << exit(FatalIOError);
        }

        const word unitSetCoeffs(word(dict.lookup("unitSet")) + "Coeffs");

        if (!dict.found(unitSetCoeffs))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find " << unitSetCoeffs << " in dictionary "
                << dict.name() << exit(FatalIOError);
        }

        const dictionary& unitDict = dict.subDict(unitSetCoeffs);

        unitSetPtr_ = new HashTable<dimensionedScalar>(unitDict.size());

        forAllConstIter(dictionary, unitDict, iter)
        {
            if (iter().keyword() != "writeUnits")
            {
                dimensionedScalar dt(iter().keyword(), iter().stream());
                const bool ok = unitSetPtr_->insert(iter().keyword(), dt);
                if (!ok)
                {
                    FatalIOErrorInFunction(dict)
                        << "Duplicate unit " << iter().keyword()
                        << " in DimensionSets dictionary"
                        << exit(FatalIOError);
                }
            }
        }

        const wordList writeUnitNames
        (
            unitDict.lookupOrDefault<wordList>
            (
                "writeUnits",
                wordList(0)
            )
        );

        writeUnitSetPtr_ = new dimensionSets(*unitSetPtr_, writeUnitNames);

        if (writeUnitNames.size() != 0 && writeUnitNames.size() != 7)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find entry \"writeUnits\" in " << unitDict.name()
                << " or it is not a wordList of size 7"
                << exit(FatalIOError);
        }
    }

    return *unitSetPtr_;
}


const Foam::dimensionSets& Foam::writeUnitSet()
{
    if (!writeUnitSetPtr_)
    {
        (void)unitSet();
    }
    return *writeUnitSetPtr_;
}


const Foam::dimensionSet Foam::dimless(0, 0, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimMass(1, 0, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimLength(0, 1, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTime(0, 0, 1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const Foam::dimensionSet Foam::dimMoles(0, 0, 0, 0, 1, 0, 0);
const Foam::dimensionSet Foam::dimCurrent(0, 0, 0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const Foam::dimensionSet Foam::dimArea(sqr(dimLength));
const Foam::dimensionSet Foam::dimVolume(pow3(dimLength));
const Foam::dimensionSet Foam::dimVol(dimVolume);

const Foam::dimensionSet Foam::dimVelocity(dimLength/dimTime);
const Foam::dimensionSet Foam::dimAcceleration(dimVelocity/dimTime);

const Foam::dimensionSet Foam::dimDensity(dimMass/dimVolume);
const Foam::dimensionSet Foam::dimForce(dimMass*dimAcceleration);
const Foam::dimensionSet Foam::dimEnergy(dimForce*dimLength);
const Foam::dimensionSet Foam::dimPower(dimEnergy/dimTime);

const Foam::dimensionSet Foam::dimPressure(dimForce/dimArea);
const Foam::dimensionSet Foam::dimCompressibility(dimDensity/dimPressure);
const Foam::dimensionSet Foam::dimGasConstant(dimEnergy/dimMass/dimTemperature);
const Foam::dimensionSet Foam::dimSpecificHeatCapacity(dimGasConstant);
const Foam::dimensionSet Foam::dimViscosity(dimArea/dimTime);
const Foam::dimensionSet Foam::dimDynamicViscosity(dimDensity*dimViscosity);

const Foam::dimensionSet Foam::dimFlux(dimArea*dimVelocity);
const Foam::dimensionSet Foam::dimMassFlux(dimDensity*dimFlux);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dimensionSets::dimensionSets
(
    const HashTable<dimensionedScalar>& units,
    const wordList& unitNames
)
:
    units_(unitNames.size()),
    conversion_(unitNames.size()),
    conversionPivots_(unitNames.size()),
    valid_(false)
{
    forAll(unitNames, i)
    {
        units_.set
        (
            i,
            new dimensionedScalar
            (
                units[unitNames[i]]
            )
        );
    }

    if (unitNames.size() == 7)
    {
        valid_ = true;

        // Determine conversion from basic units to write units
        for (label rowI = 0; rowI < conversion_.m(); rowI++)
        {
            scalar* row = conversion_[rowI];

            for (label columnI = 0; columnI < conversion_.n(); columnI++)
            {
                const dimensionedScalar& dSet = units_[columnI];
                row[columnI] = dSet.dimensions()[rowI];
            }
        }

        conversionPivots_.setSize(conversion_.m());
        LUDecompose(conversion_, conversionPivots_);
    }
}


void Foam::dimensionSets::coefficients(scalarField& exponents) const
{
    LUBacksubstitute(conversion_, conversionPivots_, exponents);
}


// ************************************************************************* //
