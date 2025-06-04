/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "populationBalanceModel.H"
#include "populationBalanceSizeDistribution.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(populationBalanceSizeDistribution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        populationBalanceSizeDistribution,
        dictionary
    );
}
}

const Foam::NamedEnum
<
    Foam::functionObjects::populationBalanceSizeDistribution::functionType,
    6
>
Foam::functionObjects::populationBalanceSizeDistribution::functionTypeNames_
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
    Foam::functionObjects::populationBalanceSizeDistribution::coordinateType,
    4
>
Foam::functionObjects::populationBalanceSizeDistribution:: coordinateTypeNames_
{
    "volume",
    "area",
    "diameter",
    "projectedAreaDiameter"
};

const Foam::NamedEnum
<
    Foam::functionObjects::populationBalanceSizeDistribution::weightType,
    4
>
Foam::functionObjects::populationBalanceSizeDistribution::weightTypeNames_
{
    "numberConcentration",
    "volumeConcentration",
    "areaConcentration",
    "cellVolume"
};

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::word
Foam::functionObjects::populationBalanceSizeDistribution::
functionTypeSymbolicName()
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


Foam::word
Foam::functionObjects::populationBalanceSizeDistribution::
coordinateTypeSymbolicName
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
Foam::functionObjects::populationBalanceSizeDistribution::filterField
(
    const scalarField& field
) const
{
    if (zone_.all())
    {
        return field;
    }
    else
    {
        return tmp<scalarField>(new scalarField(field, zone_.zone()));
    }
}


Foam::scalar
Foam::functionObjects::populationBalanceSizeDistribution::averageCoordinateValue
(
    const populationBalanceModel& popBal,
    const label i,
    const coordinateType& cType
)
{
    scalar averageCoordinateValue(Zero);

    switch (cType)
    {
        case coordinateType::volume:
        {
            averageCoordinateValue = popBal.v(i).value();

            break;
        }
        case coordinateType::area:
        {
            averageCoordinateValue = weightedAverage(popBal, i, popBal.a(i));

            break;
        }
        case coordinateType::diameter:
        {
            averageCoordinateValue = weightedAverage(popBal, i, popBal.d(i));

            break;
        }
        case coordinateType::projectedAreaDiameter:
        {
            averageCoordinateValue =
                weightedAverage(popBal, i, sqrt(popBal.a(i)/pi));

            break;
        }
    }

    return averageCoordinateValue;
}


Foam::scalar
Foam::functionObjects::populationBalanceSizeDistribution::weightedAverage
(
    const populationBalanceModel& popBal,
    const label i,
    const Foam::scalarField& field
)
{
    const volScalarField& alpha = popBal.phases()[i];
    const volScalarField& fi = popBal.f(i);

    scalar weightedAverage(Zero);

    switch (weightType_)
    {
        case weightType::numberConcentration:
        {
            scalarField Ni(filterField(fi*alpha/popBal.v(i).value()));

            if (gSum(Ni) == 0)
            {
                weightedAverage =
                    gSum(filterField(mesh_.V()*field))/zone_.V();
            }
            else
            {
                weightedAverage =
                    gSum(Ni*filterField(field))/gSum(Ni);
            }

            break;
        }
        case weightType::volumeConcentration:
        {
            scalarField Vi(filterField(fi*alpha));

            if (gSum(Vi) == 0)
            {
                weightedAverage =
                    gSum(filterField(mesh_.V()*field))/zone_.V();
            }
            else
            {
                weightedAverage =
                    gSum(Vi*filterField(field))/gSum(Vi);
            }

            break;
        }
        case weightType::areaConcentration:
        {
            scalarField Ai(filterField(popBal.a(i)()*alpha));

            if (gSum(Ai) == 0)
            {
                weightedAverage =
                    gSum(filterField(mesh_.V()*field))/zone_.V();
            }
            else
            {
                weightedAverage =
                    gSum(Ai*filterField(field))/gSum(Ai);
            }

            break;
        }
        case weightType::cellVolume:
        {
            weightedAverage =
                gSum(filterField(mesh_.V()*field))/zone_.V();

            break;
        }
    }

    return weightedAverage;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalanceSizeDistribution::
populationBalanceSizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    file_(obr_, name),
    mesh_(fvMeshFunctionObject::mesh_),
    zone_(fvMeshFunctionObject::mesh_, dict),
    popBalName_(dict.lookup("populationBalance")),
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

Foam::functionObjects::populationBalanceSizeDistribution::
~populationBalanceSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::populationBalanceSizeDistribution::read
(
    const dictionary& dict
)
{
    Log << type() << " " << name() << ":" << nl;

    fvMeshFunctionObject::read(dict);

    formatterPtr_ = setWriter::New(dict.lookup("setFormat"), dict);

    return false;
}


bool Foam::functionObjects::populationBalanceSizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::populationBalanceSizeDistribution::write()
{
    Log << type() << " " << name() << " write:" << nl;

    const populationBalanceModel& popBal =
        obr_.lookupObject<populationBalanceModel>(popBalName_);

    scalarField coordinateValues(popBal.nGroups());
    scalarField boundaryValues(popBal.nGroups() + 1);
    scalarField resultValues(popBal.nGroups());

    forAll(popBal.fs(), i)
    {
        coordinateValues[i] =
            averageCoordinateValue(popBal, i, coordinateType_);
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
            forAll(popBal.fs(), i)
            {
                const volScalarField& alpha = popBal.phases()[i];
                const volScalarField& fi = popBal.f(i);
                const dimensionedScalar& vi = popBal.v(i);

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*alpha/vi))/zone_.V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            break;
        }
        case functionType::numberDensity:
        {
            forAll(popBal.fs(), i)
            {
                const volScalarField& alpha = popBal.phases()[i];
                const volScalarField& fi = popBal.f(i);
                const dimensionedScalar& vi = popBal.v(i);

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*alpha/vi))/zone_.V();
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
            forAll(popBal.fs(), i)
            {
                const volScalarField& alpha = popBal.phases()[i];
                const volScalarField& fi = popBal.f(i);

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*alpha))/zone_.V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            break;
        }
        case functionType::volumeDensity:
        {
            forAll(popBal.fs(), i)
            {
                const volScalarField& alpha = popBal.phases()[i];
                const volScalarField& fi = popBal.f(i);

                resultValues[i] =
                    gSum(filterField(mesh_.V()*fi*alpha))/zone_.V();
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
            forAll(popBal.fs(), i)
            {
                const volScalarField& alpha = popBal.phases()[i];
                const volScalarField& fi = popBal.f(i);
                const dimensionedScalar& vi = popBal.v(i);

                resultValues[i] =
                    gSum
                    (
                        filterField(mesh_.V()*popBal.a(i)()*fi*alpha/vi)
                    )
                   /zone_.V();
            }

            if (normalise_ && sum(resultValues) != 0)
            {
                resultValues /= sum(resultValues);
            }

            break;
        }
        case functionType::areaDensity:
        {
            forAll(popBal.fs(), i)
            {
                const volScalarField& alpha = popBal.phases()[i];
                const volScalarField& fi = popBal.f(i);
                const dimensionedScalar& vi = popBal.v(i);

                resultValues[i] =
                    gSum
                    (
                        filterField(mesh_.V()*popBal.a(i)()*fi*alpha/vi)
                    )
                   /zone_.V();
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

        forAll(coordinateTypeNames_, i)
        {
            const coordinateType cType = coordinateType(i);

            otherCoordinateSymbolicNames[cType] =
                coordinateTypeSymbolicName(cType);

            otherCoordinateValues.set(cType, new scalarField(popBal.nGroups()));

            forAll(popBal.fs(), i)
            {
                otherCoordinateValues[cType][i] =
                    averageCoordinateValue(popBal, i, cType);
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


void Foam::functionObjects::populationBalanceSizeDistribution::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &this->mesh())
    {
        zone_.movePoints();
    }
}


void Foam::functionObjects::populationBalanceSizeDistribution::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.topoChange(map);
    }
}


void Foam::functionObjects::populationBalanceSizeDistribution::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.mapMesh(map);
    }
}


void Foam::functionObjects::populationBalanceSizeDistribution::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.distribute(map);
    }
}


// ************************************************************************* //
