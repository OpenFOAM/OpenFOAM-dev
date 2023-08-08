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
const Foam::dimensionSet Foam::dimMomentum(dimMass*dimVelocity);
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

dictionary* dimensionSetsDictPtr_(nullptr);

const dictionary& dimensionSetsDict()
{
    if (!dimensionSetsDictPtr_)
    {
        dictionary* cachedPtr = nullptr;
        dimensionSetsDictPtr_ = new dictionary
        (
            debug::switchSet
            (
                "DimensionSets",
                cachedPtr
            )
        );
    }

    return *dimensionSetsDictPtr_;
}

HashTable<dimensionedScalar>* addedUnitsPtr_(nullptr);
HashTable<dimensionedScalar>* unitSetPtr_(nullptr);

// Delete the above data at the end of the run
struct deleteDimensionSystemsPtr
{
    ~deleteDimensionSystemsPtr()
    {
        deleteDemandDrivenData(dimensionSetsDictPtr_);
        deleteDemandDrivenData(addedUnitsPtr_);
        deleteDemandDrivenData(unitSetPtr_);
    }
};

deleteDimensionSystemsPtr deleteDimensionSystemsPtr_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::addUnit(const dimensionedScalar& unit)
{
    deleteDemandDrivenData(dimensionSetsDictPtr_);

    if (!addedUnitsPtr_)
    {
        addedUnitsPtr_ = new HashTable<dimensionedScalar>();
    }

    addedUnitsPtr_->insert(unit.name(), unit);

    deleteDemandDrivenData(unitSetPtr_);
}


const Foam::HashTable<Foam::dimensionedScalar>& Foam::unitSet()
{
    if (!unitSetPtr_)
    {
        const dictionary& dimSetsDict = dimensionSetsDict();

        if (!dimSetsDict.found("unitSet"))
        {
            FatalIOErrorInFunction(dimSetsDict)
                << "Cannot find unitSet in dictionary " << dimSetsDict.name()
                << exit(FatalIOError);
        }

        const word unitSetDictName =
            dimSetsDict.lookup<word>("unitSet") + "Coeffs";

        if (!dimSetsDict.found(unitSetDictName))
        {
            FatalIOErrorInFunction(dimSetsDict)
                << "Cannot find " << unitSetDictName << " in dictionary "
                << dimSetsDict.name() << exit(FatalIOError);
        }

        const dictionary& unitSetDict = dimSetsDict.subDict(unitSetDictName);

        unitSetPtr_ = new HashTable<dimensionedScalar>(unitSetDict.size());

        forAllConstIter(dictionary, unitSetDict, iter)
        {
            const dimensionedScalar dt(iter().keyword(), iter().stream());

            const bool ok = unitSetPtr_->insert(iter().keyword(), dt);

            if (!ok)
            {
                FatalIOErrorInFunction(dimSetsDict)
                    << "Duplicate unit " << iter().keyword()
                    << " read from dictionary"
                    << exit(FatalIOError);
            }
        }

        if (addedUnitsPtr_)
        {
            forAllConstIter(HashTable<dimensionedScalar>, *addedUnitsPtr_, iter)
            {
                const bool ok = unitSetPtr_->insert(iter.key(), iter());

                if (!ok)
                {
                    FatalIOErrorInFunction(dimSetsDict)
                        << "Duplicate unit " << iter.key()
                        << " added to dictionary"
                        << exit(FatalIOError);
                }
            }
        }
    }

    return *unitSetPtr_;
}


// ************************************************************************* //
