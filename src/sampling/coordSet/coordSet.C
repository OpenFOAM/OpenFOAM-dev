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

#include "coordSet.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum<Foam::coordSet::axisType, 6>::names[] =
    {
        "xyz",
        "x",
        "y",
        "z",
        "distance",
        "default"
    };
}


const Foam::NamedEnum<Foam::coordSet::axisType, 6>
    Foam::coordSet::axisTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSet::coordSet()
:
    segments_(0),
    positionName_(word::null),
    positions_(nullptr),
    distanceName_(word::null),
    distances_(nullptr),
    axis_(axisType::DEFAULT)
{}


Foam::coordSet::coordSet
(
    const labelList& segments,
    const word& positionName,
    const pointField& positions,
    const word& distanceName,
    const scalarField& distances,
    const word& axis
)
:
    segments_(segments),
    positionName_(positionName),
    positions_
    (
        isNull<pointField>(positions)
      ? nullptr
      : new pointField(positions)
    ),
    distanceName_(distanceName),
    distances_
    (
        isNull<scalarField>(distances)
      ? nullptr
      : new scalarField(distances)
    ),
    axis_(axisTypeNames_[axis])
{}


Foam::coordSet::coordSet
(
    const bool contiguous,
    const word& positionName,
    const pointField& positions,
    const word& axis
)
:
    coordSet
    (
        contiguous
      ? labelList(positions.size(), 0)
      : identityMap(positions.size()),
        positionName,
        positions,
        word::null,
        scalarField::null(),
        axis
    )
{}


Foam::coordSet::coordSet
(
    const bool contiguous,
    const word& distanceName,
    const scalarField& distances,
    const word& axis
)
:
    coordSet
    (
        contiguous
      ? labelList(distances.size(), 0)
      : identityMap(distances.size()),
        word::null,
        pointField::null(),
        distanceName,
        distances,
        axis
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSet::hasScalarAxis() const
{
    // If axis is default and both positions and distances are valid, take
    // scalar distances in preference

    return
        (axis_ == axisType::X && positions_.valid())
     || (axis_ == axisType::Y && positions_.valid())
     || (axis_ == axisType::Z && positions_.valid())
     || (axis_ == axisType::DISTANCE && distances_.valid())
     || (axis_ == axisType::DEFAULT && distances_.valid());
}


bool Foam::coordSet::hasPointAxis() const
{
    // If axis is default and both positions and distances are valid, take
    // scalar distances in preference

    return
        (axis_ == axisType::XYZ && positions_.valid())
     || (axis_ == axisType::DEFAULT && positions_.valid());
}


Foam::scalar Foam::coordSet::scalarCoord(const label index) const
{
    switch (axis_)
    {
        case axisType::XYZ:
            FatalErrorInFunction
                << "Scalar coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis" << exit(FatalError);
            break;
        case axisType::X:
            return positions_()[index].x();
        case axisType::Y:
            return positions_()[index].y();
        case axisType::Z:
            return positions_()[index].z();
        case axisType::DISTANCE:
            return distances_()[index];
        case axisType::DEFAULT:
            if (distances_.valid()) return distances_()[index];
            FatalErrorInFunction
                << "Scalar coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis, but with no valid distances"
                << exit(FatalError);
            break;
    }

    return NaN;
}


Foam::tmp<Foam::scalarField> Foam::coordSet::scalarCoords() const
{
    switch (axis_)
    {
        case axisType::XYZ:
            FatalErrorInFunction
                << "Scalar coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis" << exit(FatalError);
            break;
        case axisType::X:
            return positions_().component(point::X);
        case axisType::Y:
            return positions_().component(point::Y);
        case axisType::Z:
            return positions_().component(point::Z);
        case axisType::DISTANCE:
            return distances_();
        case axisType::DEFAULT:
            if (distances_.valid()) return distances_();
            FatalErrorInFunction
                << "Scalar coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis, but with no valid distances"
                << exit(FatalError);
            break;
    }

    return scalarField::null();
}


Foam::word Foam::coordSet::scalarName() const
{
    switch (axis_)
    {
        case axisType::XYZ:
            FatalErrorInFunction
                << "Scalar name requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis" << exit(FatalError);
            break;
        case axisType::X:
        case axisType::Y:
        case axisType::Z:
            return
                positionName_
              + (positionName_ != word::null ? "_" : "")
              + axis();
        case axisType::DISTANCE:
            return distanceName_;
        case axisType::DEFAULT:
            if (distances_.valid()) return distanceName_;
            FatalErrorInFunction
                << "Scalar coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis, but with no valid distances"
                << exit(FatalError);
            break;
    }

    return word::null;
}


Foam::point Foam::coordSet::pointCoord(const label index) const
{
    switch (axis_)
    {
        case axisType::XYZ:
            return positions_()[index];
        case axisType::X:
        case axisType::Y:
        case axisType::Z:
        case axisType::DISTANCE:
            FatalErrorInFunction
                << "Point coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis" << exit(FatalError);
            break;
        case axisType::DEFAULT:
            if (positions_.valid()) return positions_()[index];
            FatalErrorInFunction
                << "Point coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis, but with no valid point"
                << exit(FatalError);
            break;
    }

    return point::uniform(NaN);
}


Foam::tmp<Foam::pointField> Foam::coordSet::pointCoords() const
{
    switch (axis_)
    {
        case axisType::XYZ:
            return positions_();
        case axisType::X:
        case axisType::Y:
        case axisType::Z:
        case axisType::DISTANCE:
            FatalErrorInFunction
                << "Point coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis" << exit(FatalError);
            break;
        case axisType::DEFAULT:
            if (positions_.valid()) return positions_();
            FatalErrorInFunction
                << "Point coordinate requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis, but with no valid point"
                << exit(FatalError);
            break;
    }

    return pointField::null();
}


Foam::word Foam::coordSet::pointName() const
{
    switch (axis_)
    {
        case axisType::XYZ:
            return positionName_;
        case axisType::X:
        case axisType::Y:
        case axisType::Z:
        case axisType::DISTANCE:
            FatalErrorInFunction
                << "Point name requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis" << exit(FatalError);
            break;
        case axisType::DEFAULT:
            if (positions_.valid()) return positionName_;
            FatalErrorInFunction
                << "Point name requested from coordinate set with "
                << axisTypeNames_[axis_] << " axis, but with no valid point"
                << exit(FatalError);
            break;
    }

    return word::null;
}


Foam::labelList Foam::coordSet::vertices() const
{
    label nVertices = 0;
    labelList vertices(size());

    forAll(*this, pointi)
    {
        const label s = segments_[pointi];
        if
        (
            (pointi == 0 || s != segments_[pointi - 1])
         && (pointi == size() - 1 || s != segments_[pointi + 1])
        )
        {
            vertices[nVertices ++] = pointi;
        }
    }

    vertices.resize(nVertices);

    return vertices;
}


Foam::labelPairList Foam::coordSet::edges() const
{
    label nEdges = 0;
    labelPairList edges(size() - 1);

    for (label pointi = 0; pointi < size() - 1; ++ pointi)
    {
        if (segments_[pointi] == segments_[pointi + 1])
        {
            edges[nEdges ++] = labelPair(pointi, pointi + 1);
        }
    }

    edges.resize(nEdges);

    return edges;
}


Foam::labelListList Foam::coordSet::lines() const
{
    DynamicList<label> line;
    DynamicList<labelList> lines;

    forAll(*this, pointi)
    {
        line.append(pointi);

        if
        (
            pointi == size() - 1
         || segments_[pointi] != segments_[pointi + 1]
        )
        {
            if (line.size() > 1)
            {
                lines.append(line);
            }

            line.clear();
        }
    }

    labelListList linesNonDynamic;
    linesNonDynamic.transfer(lines);

    return linesNonDynamic;
}


Foam::Tuple2<Foam::coordSet, Foam::labelList> Foam::coordSet::gather() const
{
    // Collect data from all processors
    List<List<point>> gatheredPositions;
    if (positions_.valid())
    {
        gatheredPositions.resize(Pstream::nProcs());
        gatheredPositions[Pstream::myProcNo()] = positions_();
        Pstream::gatherList(gatheredPositions);
    }

    List<scalarList> gatheredDistances;
    if (distances_.valid())
    {
        gatheredDistances.resize(Pstream::nProcs());
        gatheredDistances[Pstream::myProcNo()] = distances_();
        Pstream::gatherList(gatheredDistances);
    }

    List<labelField> gatheredSegments(Pstream::nProcs());
    gatheredSegments[Pstream::myProcNo()] = segments_;
    Pstream::gatherList(gatheredSegments);

    // Combine processor lists into one big list.
    List<point> allPositions;
    if (positions_.valid())
    {
        allPositions =
            ListListOps::combine<List<point>>
            (
                gatheredPositions,
                accessOp<List<point>>()
            );
    };

    scalarList allDistances;
    if (distances_.valid())
    {
        allDistances =
            ListListOps::combine<scalarList>
            (
                gatheredDistances,
                accessOp<scalarList>()
            );
    }

    labelList allSegments
    (
        ListListOps::combine<labelList>
        (
            gatheredSegments,
            accessOp<labelList>()
        )
    );

    // Construct a result tuple
    Tuple2<coordSet, labelList> result
    (
        coordSet(),
        identityMap(allSegments.size())
    );
    coordSet& set = result.first();
    labelList& order = result.second();

    // Sort by segment then by distance
    if (distances_.valid())
    {
        stableSort
        (
            order,
            [&](const label a, const label b)
            {
                return
                    allSegments[a] < allSegments[b] ? true
                  : allSegments[a] > allSegments[b] ? false
                  : allDistances[a] < allDistances[b];
            }
        );
    }
    else
    {
        stableSort
        (
            order,
            [&](const label a, const label b)
            {
                return allSegments[a] < allSegments[b];
            }
        );
    }

    // Set the data in the coordinate set
    set.segments_ = labelField(allSegments, order);
    set.positionName_ = positionName_;
    if (positions_.valid())
    {
        set.positions_.set(new pointField(allPositions, order));
    }
    set.distanceName_ = distanceName_;
    if (distances_.valid())
    {
        set.distances_.set(new scalarField(allDistances, order));
    }
    set.axis_ = axis_;

    return result;
}


// ************************************************************************* //
