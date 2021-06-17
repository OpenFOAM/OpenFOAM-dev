/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "meshSearch.H"
#include "setWriter.H"
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

void Foam::sampledSet::setSamples
(
    const List<point>& samplingPts,
    const labelList& samplingCells,
    const labelList& samplingFaces,
    const labelList& samplingSegments,
    const scalarList& samplingCurveDist
)
{
    setSize(samplingPts.size());
    cells_.setSize(samplingCells.size());
    faces_.setSize(samplingFaces.size());
    segments_.setSize(samplingSegments.size());
    curveDist_.setSize(samplingCurveDist.size());

    if
    (
        (cells_.size() != size())
     || (faces_.size() != size())
     || (segments_.size() != size())
     || (curveDist_.size() != size())
    )
    {
        FatalErrorInFunction
            << "sizes not equal : "
            << "  points:" << size()
            << "  cells:" << cells_.size()
            << "  faces:" << faces_.size()
            << "  segments:" << segments_.size()
            << "  curveDist:" << curveDist_.size()
            << abort(FatalError);
    }

    forAll(samplingPts, sampleI)
    {
        operator[](sampleI) = samplingPts[sampleI];
    }
    curveDist_ = samplingCurveDist;

    cells_ = samplingCells;
    faces_ = samplingFaces;
    segments_ = samplingSegments;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSet::sampledSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis
)
:
    coordSet(name, axis),
    mesh_(mesh),
    searchEngine_(searchEngine),
    segments_(0),
    cells_(0),
    faces_(0)
{}


Foam::sampledSet::sampledSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    coordSet(name, dict.lookup("axis")),
    mesh_(mesh),
    searchEngine_(searchEngine),
    segments_(0),
    cells_(0),
    faces_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSet::~sampledSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sampledSet> Foam::sampledSet::New
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
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
            searchEngine,
            dict.optionalSubDict(sampleType + "Coeffs")
        )
    );
}


Foam::Ostream& Foam::sampledSet::write(Ostream& os) const
{
    coordSet::write(os);

    os  << endl << "\t(celli)\t(facei)" << endl;

    forAll(*this, sampleI)
    {
        os  << '\t' << cells_[sampleI]
            << '\t' << faces_[sampleI]
            << endl;
    }

    return os;
}


// ************************************************************************* //
