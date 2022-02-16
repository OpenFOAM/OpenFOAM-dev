/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

#include "sizeDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sizeDistribution, 0);
    addToRunTimeSelectionTable(functionObject, sizeDistribution, dictionary);
}
}

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::functionType,
    6
>::names[] =
{
    "numberConcentration",
    "numberDensity",
    "volumeConcentration",
    "volumeDensity",
    "areaConcentration",
    "areaDensity"
};

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::functionType,
    6
> Foam::functionObjects::sizeDistribution::functionTypeNames_;

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::coordinateType,
    4
>::names[] =
{
    "volume",
    "area",
    "diameter",
    "projectedAreaDiameter"
};

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::coordinateType,
    4
> Foam::functionObjects::sizeDistribution::coordinateTypeNames_;


namespace Foam
{
    template<>
    const char* NamedEnum
    <
        Foam::functionObjects::sizeDistribution::weightType,
        4
    >::names[] =
    {
        "numberConcentration",
        "volumeConcentration",
        "areaConcentration",
        "cellVolume"
    };
}


const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::weightType,
    4
>
Foam::functionObjects::sizeDistribution::weightTypeNames_;

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::word Foam::functionObjects::sizeDistribution::functionTypeSymbolicName()
{
    word functionTypeSymbolicName(word::null);

    switch (functionType_)
    {
        case functionType::numberConcentration:
        {
            functionTypeSymbolicName = "N";

            break;
        }
        case functionType::numberDensity:
        {
            functionTypeSymbolicName = "n";

            break;
        }
        case functionType::volumeConcentration:
        {
            functionTypeSymbolicName = "V";

            break;
        }
        case functionType::volumeDensity:
        {
            functionTypeSymbolicName = "v";

            break;
        }
        case functionType::areaConcentration:
        {
            functionTypeSymbolicName = "A";

            break;
        }
        case functionType::areaDensity:
        {
            functionTypeSymbolicName = "a";

            break;
        }
    }

    return functionTypeSymbolicName;
}


Foam::word Foam::functionObjects::sizeDistribution::coordinateTypeSymbolicName
(
    const coordinateType& cType
)
{
    word coordinateTypeSymbolicName(word::null);

    switch (cType)
    {
        case coordinateType::volume:
        {
            coordinateTypeSymbolicName = "v";

            break;
        }
        case coordinateType::area:
        {
            coordinateTypeSymbolicName = "a";

            break;
        }
        case coordinateType::diameter:
        {
            coordinateTypeSymbolicName = "d";

            break;
        }
        case coordinateType::projectedAreaDiameter:
        {
            coordinateTypeSymbolicName = "dPa";

            break;
        }
    }

    return coordinateTypeSymbolicName;
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::sizeDistribution::filterField
(
    const scalarField& field
) const
{
    if (isNull(cellIDs()))
    {
        return field;
    }
    else
    {
        return tmp<scalarField>(new scalarField(field, cellIDs()));
    }
}


Foam::scalar Foam::functionObjects::sizeDistribution::averageCoordinateValue
(
    const Foam::diameterModels::sizeGroup& fi,
    const coordinateType& cType
)
{
    scalar averageCoordinateValue(Zero);

    switch (cType)
    {
        case coordinateType::volume:
        {
            averageCoordinateValue = fi.x().value();

            break;
        }
        case coordinateType::area:
        {
            averageCoordinateValue =
                weightedAverage(fi.a(), fi);

            break;
        }
        case coordinateType::diameter:
        {
            averageCoordinateValue =
                weightedAverage(fi.d(), fi);

            break;
        }
        case coordinateType::projectedAreaDiameter:
        {
            averageCoordinateValue =
                weightedAverage(sqrt(fi.a()/pi), fi);

            break;
        }
    }

    return averageCoordinateValue;
}


Foam::scalar Foam::functionObjects::sizeDistribution::weightedAverage
(
    const Foam::scalarField& fld,
    const Foam::diameterModels::sizeGroup& fi
)
{
    scalar weightedAverage(Zero);

    switch (weightType_)
    {
        case weightType::numberConcentration:
        {
            scalarField Ni(filterField(fi*fi.phase()/fi.x().value()));

            if (gSum(Ni) == 0)
            {
                weightedAverage =
                    gSum(filterField(mesh_.V()*fld))/this->V();
            }
            else
            {
                weightedAverage =
                    gSum(Ni*filterField(fld))/gSum(Ni);
            }

            break;
        }
        case weightType::volumeConcentration:
        {
            scalarField Vi(filterField(fi*fi.phase()));

            if (gSum(Vi) == 0)
            {
                weightedAverage =
                    gSum(filterField(mesh_.V()*fld))/this->V();
            }
            else
            {
                weightedAverage =
                    gSum(Vi*filterField(fld))/gSum(Vi);
            }

            break;
        }
        case weightType::areaConcentration:
        {
            scalarField Ai(filterField(fi.a().ref()*fi.phase()));

            if (gSum(Ai) == 0)
            {
                weightedAverage =
                    gSum(filterField(mesh_.V()*fld))/this->V();
            }
            else
            {
                weightedAverage =
                    gSum(Ai*filterField(fld))/gSum(Ai);
            }

            break;
        }
        case weightType::cellVolume:
        {
            weightedAverage =
                gSum(filterField(mesh_.V()*fld))/this->V();

            break;
        }
    }

    return weightedAverage;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::sizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    file_(obr_, name),
    mesh_(fvMeshFunctionObject::mesh_),
    popBal_
    (
        obr_.lookupObject<Foam::diameterModels::populationBalanceModel>
        (
            dict.lookup("populationBalance")
        )
    ),
    functionType_(functionTypeNames_.read(dict.lookup("functionType"))),
    coordinateType_(coordinateTypeNames_.read(dict.lookup("coordinateType"))),
    allCoordinates_
    (
        dict.lookupOrDefault<Switch>("allCoordinates", false)
    ),
    normalise_(dict.lookupOrDefault<Switch>("normalise", false)),
    logTransform_
    (
        dict.lookupOrDefaultBackwardsCompatible<Switch>
        (
            {"logTransform", "geometric"},
            false
        )
    ),
    weightType_
    (
        dict.found("weightType")
      ? weightTypeNames_.read(dict.lookup("weightType"))
      : weightType::numberConcentration
    ),
    formatterPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::~sizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sizeDistribution::read(const dictionary& dict)
{
    Log << type() << " " << name() << ":" << nl;

    fvMeshFunctionObject::read(dict);

    formatterPtr_ = setWriter::New(dict.lookup("setFormat"), dict);

    return false;
}


bool Foam::functionObjects::sizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::sizeDistribution::write()
{
    Log << type() << " " << name() << " write:" << nl;

    const UPtrList<diameterModels::sizeGroup>& sizeGroups =
        popBal_.sizeGroups();

    scalarField coordinateValues(sizeGroups.size());
    scalarField boundaryValues(sizeGroups.size() + 1);
    scalarField resultValues(sizeGroups.size());

    forAll(sizeGroups, i)
    {
        const diameterModels::sizeGroup& fi = sizeGroups[i];

        coordinateValues[i] = averageCoordinateValue(fi, coordinateType_);
    }

    if
    (
        functionType_ == functionType::numberDensity
     || functionType_ == functionType::volumeDensity
     || functionType_ == functionType::areaDensity
    )
    {
        boundaryValues.first() = coordinateValues.first();
        boundaryValues.last() = coordinateValues.last();

        for (label i = 1; i < boundaryValues.size() - 1; i++)
        {
            boundaryValues[i] =
                0.5*(coordinateValues[i] + coordinateValues[i-1]);
        }

        if (logTransform_)
        {
            boundaryValues = Foam::log(boundaryValues);
        }
    }

    switch (functionType_)
    {
        case functionType::numberConcentration:
        {
            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*fi.phase()/fi.x()))/this->V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            break;
        }
        case functionType::numberDensity:
        {
            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*fi.phase()/fi.x()))/this->V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            forAll(resultValues, i)
            {
                resultValues[i] /= (boundaryValues[i+1] - boundaryValues[i]);
            }

            break;
        }
        case functionType::volumeConcentration:
        {
            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*fi.phase()))/this->V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            break;
        }
        case functionType::volumeDensity:
        {
            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*fi.phase()))/this->V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            forAll(resultValues, i)
            {
                resultValues[i] /= (boundaryValues[i+1] - boundaryValues[i]);
            }

            break;
        }
        case functionType::areaConcentration:
        {
            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                resultValues[i] =
                    gSum
                    (
                        filterField(mesh_.V()*fi.a().ref()*fi*fi.phase()/fi.x())
                    )
                   /this->V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            break;
        }
        case functionType::areaDensity:
        {
            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                resultValues[i] =
                    gSum
                    (
                        filterField(mesh_.V()*fi.a().ref()*fi*fi.phase()/fi.x())
                    )
                   /this->V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            forAll(resultValues, i)
            {
                resultValues[i] /= (boundaryValues[i+1] - boundaryValues[i]);
            }

            break;
        }
    }


    if (allCoordinates_)
    {
        wordList otherCoordinateSymbolicNames(coordinateTypeNames_.size());
        PtrList<scalarField> otherCoordinateValues(coordinateTypeNames_.size());
        typedef NamedEnum<coordinateType, 4> namedEnumCoordinateType;

        forAllConstIter(namedEnumCoordinateType, coordinateTypeNames_, iter)
        {
            const coordinateType cType = coordinateTypeNames_[iter.key()];

            otherCoordinateSymbolicNames[cType] =
                coordinateTypeSymbolicName(cType);

            otherCoordinateValues.set
            (
                cType,
                new scalarField(popBal_.sizeGroups().size())
            );

            forAll(sizeGroups, i)
            {
                const diameterModels::sizeGroup& fi = sizeGroups[i];

                otherCoordinateValues[cType][i] =
                    averageCoordinateValue(fi, cType);
            }
        }

        if (Pstream::master())
        {
            formatterPtr_->write
            (
                file_.baseTimeDir(),
                name(),
                coordSet
                (
                    true,
                    coordinateTypeSymbolicName(coordinateType_),
                    coordinateValues
                ),
                functionTypeSymbolicName(),
                resultValues,
                otherCoordinateSymbolicNames,
                otherCoordinateValues
            );
        }
    }
    else
    {
        if (Pstream::master())
        {
            formatterPtr_->write
            (
                file_.baseTimeDir(),
                name(),
                coordSet
                (
                    true,
                    coordinateTypeSymbolicName(coordinateType_),
                    coordinateValues
                ),
                functionTypeSymbolicName(),
                resultValues
            );
        }
    }

    return true;
}


// ************************************************************************* //
