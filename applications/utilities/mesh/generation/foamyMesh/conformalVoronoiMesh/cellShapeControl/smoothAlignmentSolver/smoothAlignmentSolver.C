/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "smoothAlignmentSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Triangulation, class Type>
Foam::tmp<Foam::Field<Type>> Foam::smoothAlignmentSolver::filterFarPoints
(
    const Triangulation& mesh,
    const Field<Type>& field
)
{
    tmp<Field<Type>> tNewField(new Field<Type>(field.size()));
    Field<Type>& newField = tNewField();

    label added = 0;
    label count = 0;

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            newField[added++] = field[count];
        }

        count++;
    }

    newField.resize(added);

    return tNewField;
}


template<class Triangulation>
Foam::autoPtr<Foam::mapDistribute> Foam::smoothAlignmentSolver::buildReferredMap
(
    const Triangulation& mesh,
    labelList& indices
)
{
    globalIndex globalIndexing(mesh.vertexCount());

    DynamicList<label> dynIndices(mesh.vertexCount()/10);

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (vit->referred())
        {
            dynIndices.append
            (
                globalIndexing.toGlobal(vit->procIndex(), vit->index())
            );
        }
    }

    indices.transfer(dynIndices);

    List<Map<label>> compactMap;
    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            globalIndexing,
            indices,
            compactMap
        )
    );
}


template<class Triangulation>
Foam::autoPtr<Foam::mapDistribute> Foam::smoothAlignmentSolver::buildMap
(
    const Triangulation& mesh,
    labelListList& pointPoints
)
{
    pointPoints.setSize(mesh.vertexCount());

    globalIndex globalIndexing(mesh.vertexCount());

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        std::list<typename Triangulation::Vertex_handle> adjVerts;
        mesh.finite_adjacent_vertices(vit, std::back_inserter(adjVerts));

        DynamicList<label> indices(adjVerts.size());

        for
        (
            typename std::list<typename Triangulation::Vertex_handle>::
                const_iterator adjVertI = adjVerts.begin();
            adjVertI != adjVerts.end();
            ++adjVertI
        )
        {
            typename Triangulation::Vertex_handle vh = *adjVertI;

            if (!vh->farPoint())
            {
                indices.append
                (
                    globalIndexing.toGlobal(vh->procIndex(), vh->index())
                );
            }
        }

        pointPoints[vit->index()].transfer(indices);
    }

    List<Map<label>> compactMap;
    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            globalIndexing,
            pointPoints,
            compactMap
        )
    );
}


template<class Triangulation>
Foam::tmp<Foam::triadField> Foam::smoothAlignmentSolver::buildAlignmentField
(
    const Triangulation& mesh
)
{
    tmp<triadField> tAlignments
    (
        new triadField(mesh.vertexCount(), triad::unset)
    );
    triadField& alignments = tAlignments.ref();

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        alignments[vit->index()] = vit->alignment();
    }

    return tAlignments;
}


template<class Triangulation>
Foam::tmp<Foam::pointField> Foam::smoothAlignmentSolver::buildPointField
(
    const Triangulation& mesh
)
{
    tmp<pointField> tPoints
    (
        new pointField(mesh.vertexCount(), point(great, great, great))
    );
    pointField& points = tPoints.ref();

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        points[vit->index()] = topoint(vit->point());
    }

    return tPoints;
}


void Foam::smoothAlignmentSolver::applyBoundaryConditions
(
    const triad& fixedAlignment,
    triad& t
) const
{
    label nFixed = 0;

    forAll(fixedAlignment, dirI)
    {
        if (fixedAlignment.set(dirI))
        {
            nFixed++;
        }
    }

    if (nFixed == 1)
    {
        forAll(fixedAlignment, dirI)
        {
            if (fixedAlignment.set(dirI))
            {
                t.align(fixedAlignment[dirI]);
            }
        }
    }
    else if (nFixed == 2)
    {
        forAll(fixedAlignment, dirI)
        {
            if (fixedAlignment.set(dirI))
            {
                t[dirI] = fixedAlignment[dirI];
            }
            else
            {
                t[dirI] = triad::unset[dirI];
            }
        }

        t.orthogonalize();
    }
    else if (nFixed == 3)
    {
        forAll(fixedAlignment, dirI)
        {
            if (fixedAlignment.set(dirI))
            {
                t[dirI] = fixedAlignment[dirI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothAlignmentSolver::smoothAlignmentSolver(cellShapeControlMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smoothAlignmentSolver::~smoothAlignmentSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::smoothAlignmentSolver::smoothAlignments
(
    const label maxSmoothingIterations
)
{
    scalar minResidual = 0;

    labelListList pointPoints;
    autoPtr<mapDistribute> meshDistributor = buildMap
    (
        mesh_,
        pointPoints
    );

    triadField alignments(buildAlignmentField(mesh_));
    pointField points(buildPointField(mesh_));

    // Setup the sizes and alignments on each point
    triadField fixedAlignments(mesh_.vertexCount(), triad::unset);

    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            mesh_.finite_vertices_begin();
        vit != mesh_.finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            fixedAlignments[vit->index()] = vit->alignment();
        }
    }

    Info<< nl << "Smoothing alignments" << endl;


    for (label iter = 0; iter < maxSmoothingIterations; iter++)
    {
        Info<< "Iteration " << iter;

        meshDistributor().distribute(points);
        meshDistributor().distribute(fixedAlignments);
        meshDistributor().distribute(alignments);

        scalar residual = 0;

        triadField triadAv(alignments.size(), triad::unset);

        forAll(pointPoints, pI)
        {
            const labelList& pPoints = pointPoints[pI];

            if (pPoints.empty())
            {
                continue;
            }

            triad& newTriad = triadAv[pI];

            forAll(pPoints, adjPointi)
            {
                const label adjPointIndex = pPoints[adjPointi];

                scalar dist = mag(points[pI] - points[adjPointIndex]);

                triad tmpTriad = alignments[adjPointIndex];

                for (direction dir = 0; dir < 3; dir++)
                {
                    if (tmpTriad.set(dir))
                    {
                        tmpTriad[dir] *= 1.0/(dist + small);
                    }
                }

                newTriad += tmpTriad;
            }
        }

        // Update the alignment field
        forAll(alignments, pI)
        {
            const triad& oldTriad = alignments[pI];
            triad& newTriad = triadAv[pI];

            newTriad.normalize();
            newTriad.orthogonalize();

            // Enforce the boundary conditions
            const triad& fixedAlignment = fixedAlignments[pI];

            applyBoundaryConditions
            (
                fixedAlignment,
                newTriad
            );

            newTriad = newTriad.sortxyz();

            // Residual Calculation
            for (direction dir = 0; dir < 3; ++dir)
            {
                if
                (
                    newTriad.set(dir)
                 && oldTriad.set(dir)
                 && !fixedAlignment.set(dir)
                )
                {
                    residual += diff(oldTriad, newTriad);
                }
            }

            alignments[pI] = newTriad;
        }

        reduce(residual, sumOp<scalar>());

        Info<< ", Residual = "
            << residual
              /returnReduce(points.size(), sumOp<label>())
            << endl;

        if (iter > 0 && residual <= minResidual)
        {
            break;
        }
    }

    meshDistributor().distribute(alignments);

    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            mesh_.finite_vertices_begin();
        vit != mesh_.finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            vit->alignment() = alignments[vit->index()];
        }
    }

    labelList referredPoints;
    autoPtr<mapDistribute> referredDistributor = buildReferredMap
    (
        mesh_,
        referredPoints
    );

    alignments.setSize(mesh_.vertexCount());
    referredDistributor().distribute(alignments);

    label referredI = 0;
    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            mesh_.finite_vertices_begin();
        vit != mesh_.finite_vertices_end();
        ++vit
    )
    {
        if (vit->referred())
        {
            vit->alignment() = alignments[referredPoints[referredI++]];
        }
    }
}


// ************************************************************************* //
