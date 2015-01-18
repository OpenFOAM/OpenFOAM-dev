/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "DelaunayMesh.H"
#include "labelPair.H"
#include "PrintTable.H"
#include "pointIOField.H"
#include "scalarIOField.H"
#include "labelIOField.H"
#include "pointConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::DelaunayMesh<Triangulation>::DelaunayMesh(const Time& runTime)
:
    Triangulation(),
    vertexCount_(0),
    cellCount_(0),
    runTime_(runTime)
{}


template<class Triangulation>
Foam::DelaunayMesh<Triangulation>::DelaunayMesh
(
    const Time& runTime,
    const word& meshName
)
:
    Triangulation(),
    vertexCount_(0),
    cellCount_(0),
    runTime_(runTime)
{
    Info<< "Reading " << meshName << " from " << runTime.timeName() << endl;

    pointIOField pts
    (
        IOobject
        (
            "points",
            runTime.timeName(),
            meshName/polyMesh::meshSubDir,
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (pts.headerOk())
    {
        labelIOField types
        (
            IOobject
            (
                "types",
                runTime.timeName(),
                meshName,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

// Do not read in indices
//        labelIOField indices
//        (
//            IOobject
//            (
//                "indices",
//                runTime.timeName(),
//                meshName,
//                runTime,
//                IOobject::MUST_READ,
//                IOobject::NO_WRITE
//            )
//        );

        labelIOField processorIndices
        (
            IOobject
            (
                "processorIndices",
                runTime.timeName(),
                meshName,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        List<Vb> pointsToInsert(pts.size());

        forAll(pointsToInsert, pI)
        {
            pointsToInsert[pI] =
                Vb
                (
                    toPoint(pts[pI]),
                    pI,
                    static_cast<indexedVertexEnum::vertexType>(types[pI]),
                    processorIndices[pI]
                );
        }

        rangeInsertWithInfo
        (
            pointsToInsert.begin(),
            pointsToInsert.end(),
            false,
            false
        );

        vertexCount_ = Triangulation::number_of_vertices();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::DelaunayMesh<Triangulation>::~DelaunayMesh()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::reset()
{
    Info<< "Clearing triangulation" << endl;

    DynamicList<Vb> vertices;

    for
    (
        Finite_vertices_iterator vit = Triangulation::finite_vertices_begin();
        vit != Triangulation::finite_vertices_end();
        ++vit
    )
    {
        if (vit->fixed())
        {
            vertices.append
            (
                Vb
                (
                    vit->point(),
                    vit->index(),
                    vit->type(),
                    vit->procIndex()
                )
            );

            vertices.last().fixed() = vit->fixed();
        }
    }

    this->clear();

    resetVertexCount();
    resetCellCount();

    insertPoints(vertices, false);

    Info<< "Inserted " << vertexCount() << " fixed points" << endl;
}


template<class Triangulation>
Foam::Map<Foam::label> Foam::DelaunayMesh<Triangulation>::insertPoints
(
    const List<Vb>& vertices,
    const bool reIndex
)
{
    return rangeInsertWithInfo
    (
        vertices.begin(),
        vertices.end(),
        false,
        reIndex
    );
}


template<class Triangulation>
bool Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_x_3::
operator()
(
    const Point_3& p,
    const Point_3& q
) const
{
    return typename Gt::Less_x_3()(*(p.first), *(q.first));
}

template<class Triangulation>
bool Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_y_3::
operator()
(
    const Point_3& p,
    const Point_3& q
) const
{
    return typename Gt::Less_y_3()(*(p.first), *(q.first));
}

template<class Triangulation>
bool Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_z_3::
operator()
(
    const Point_3& p,
    const Point_3& q
) const
{
    return typename Gt::Less_z_3()(*(p.first), *(q.first));
}

template<class Triangulation>
typename Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_x_3
Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::less_x_3_object()
const
{
    return Less_x_3();
}

template<class Triangulation>
typename Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_y_3
Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::less_y_3_object()
const
{
    return Less_y_3();
}

template<class Triangulation>
typename Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_z_3
Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::less_z_3_object()
const
{
    return Less_z_3();
}


template<class Triangulation>
template<class PointIterator>
Foam::Map<Foam::label> Foam::DelaunayMesh<Triangulation>::rangeInsertWithInfo
(
    PointIterator begin,
    PointIterator end,
    bool printErrors,
    bool reIndex
)
{
    typedef DynamicList
    <
        std::pair
        <
            const typename Triangulation::Point*,
            label
        >
    > vectorPairPointIndex;

    vectorPairPointIndex points;

    label count = 0;
    for (PointIterator it = begin; it != end; ++it)
    {
        points.append
        (
            std::make_pair(&(it->point()), count++)
        );
    }

    std::random_shuffle(points.begin(), points.end());

    spatial_sort
    (
        points.begin(),
        points.end(),
        Traits_for_spatial_sort()
    );

    Vertex_handle hint;

    Map<label> oldToNewIndex(points.size());

    for
    (
        typename vectorPairPointIndex::const_iterator p = points.begin();
        p != points.end();
        ++p
    )
    {
        const size_t checkInsertion = Triangulation::number_of_vertices();

        hint = this->insert(*(p->first), hint);

        const Vb& vert = *(begin + p->second);

        if (checkInsertion != Triangulation::number_of_vertices() - 1)
        {
            if (printErrors)
            {
                Vertex_handle nearV =
                    Triangulation::nearest_vertex(*(p->first));

                Pout<< "Failed insertion : " << vert.info()
                    << "         nearest : " << nearV->info();
            }
        }
        else
        {
            const label oldIndex = vert.index();
            hint->index() = getNewVertexIndex();

            if (reIndex)
            {
                oldToNewIndex.insert(oldIndex, hint->index());
            }

            hint->type() = vert.type();
            hint->procIndex() = vert.procIndex();
            hint->targetCellSize() = vert.targetCellSize();
            hint->alignment() = vert.alignment();
        }
    }

    return oldToNewIndex;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DelaunayMeshIO.C"

// ************************************************************************* //
