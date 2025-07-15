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

#include "sampledSet.H"
#include "polyMesh.H"
#include "lineCell.H"
#include "lineCellFace.H"
#include "lineFace.H"
#include "lineUniform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSet, 0);
    defineRunTimeSelectionTable(sampledSet, word);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sampledSet::setSamples() const
{
    DynamicList<point> samplingPositions;
    DynamicList<scalar> samplingDistances;
    DynamicList<label> samplingSegments;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;

    calcSamples
    (
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );

    if
    (
        (
            samplingDistances.size() != 0
         && samplingDistances.size() != samplingPositions.size()
        )
     || (samplingCells.size() != samplingPositions.size())
     || (samplingFaces.size() != samplingPositions.size())
     || (samplingSegments.size() != samplingPositions.size())
    )
    {
        FatalErrorInFunction
            << "sizes not equal : "
            << "  positions:" << samplingPositions.size()
            << "  distances:" << samplingDistances.size()
            << "  segments:" << samplingSegments.size()
            << "  cells:" << samplingCells.size()
            << "  faces:" << samplingFaces.size()
            << abort(FatalError);
    }

    pointField positions;
    positions.transfer(samplingPositions);

    scalarField distances;
    samplingDistances.transfer(samplingDistances);

    coordsPtr_.reset
    (
        new coordSet
        (
            samplingSegments,
            word::null,
            positions,
            coordSet::axisTypeNames_[coordSet::axisType::DISTANCE],
            distances.size() == positions.size()
          ? distances
          : NullObjectRef<scalarField>(),
            coordSet::axisTypeNames_[axis_]
        )
    );

    cellsPtr_.reset(new labelList());
    cellsPtr_().transfer(samplingCells);

    facesPtr_.reset(new labelList());
    facesPtr_().transfer(samplingFaces);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSet::sampledSet
(
    const word& name,
    const polyMesh& mesh,
    const word& axis
)
:
    name_(name),
    mesh_(mesh),
    coordsPtr_(nullptr),
    cellsPtr_(),
    facesPtr_(),
    axis_(coordSet::axisTypeNames_[axis])
{}


Foam::sampledSet::sampledSet
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    coordsPtr_(),
    cellsPtr_(),
    facesPtr_(),
    axis_
    (
        coordSet::axisTypeNames_
        [
            dict.lookupOrDefault<word>
            (
                "axis",
                coordSet::axisTypeNames_[coordSet::axisType::DEFAULT]
            )
        ]
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSet::~sampledSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sampledSet> Foam::sampledSet::New
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const word sampleType(dict.lookup("type"));

    const HashTable<word> oldToNewType =
    {
        Tuple2<word, word>("midPoint", sampledSets::lineCell::typeName),
        Tuple2<word, word>
        (
            "midPointAndFace",
            sampledSets::lineCellFace::typeName
        ),
        Tuple2<word, word>("face", sampledSets::lineFace::typeName),
        Tuple2<word, word>("uniform", sampledSets::lineUniform::typeName)
    };

    if (oldToNewType.found(sampleType))
    {
        const word newSampleType = oldToNewType[sampleType];

        FatalErrorInFunction
            << "Unknown sample set type "
            << sampleType << nl << nl
            << "The sample set type " << sampleType << " has been renamed "
            << newSampleType << nl << nl
            << "Replace \"type " << sampleType << ";\" with \"type "
            << newSampleType << ";\" for the set " << name << " in "
            << dict.name() << exit(FatalError);
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(sampleType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown sample set type "
            << sampleType << nl << nl
            << "Valid sample set types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<sampledSet>
    (
        cstrIter()
        (
            name,
            mesh,
            dict.optionalSubDict(sampleType + "Coeffs")
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSet::movePoints()
{
    coordsPtr_.clear();
    cellsPtr_.clear();
    facesPtr_.clear();
}


void Foam::sampledSet::topoChange(const polyTopoChangeMap& map)
{
    coordsPtr_.clear();
    cellsPtr_.clear();
    facesPtr_.clear();
}


void Foam::sampledSet::mapMesh(const polyMeshMap& map)
{
    coordsPtr_.clear();
    cellsPtr_.clear();
    facesPtr_.clear();
}


void Foam::sampledSet::distribute
(
    const polyDistributionMap& map
)
{
    coordsPtr_.clear();
    cellsPtr_.clear();
    facesPtr_.clear();
}


// ************************************************************************* //
