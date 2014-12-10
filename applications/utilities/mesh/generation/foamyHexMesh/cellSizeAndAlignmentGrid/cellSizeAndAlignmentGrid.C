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

Application
    Test-distributedDelaunayMesh

Description

\*---------------------------------------------------------------------------*/

#include "CGALTriangulation3DKernel.H"

#include "indexedVertex.H"
#include "indexedCell.H"

#include "argList.H"
#include "Time.H"
#include "DistributedDelaunayMesh.H"
#include "backgroundMeshDecomposition.H"
#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "PrintTable.H"
#include "Random.H"
#include "boundBox.H"
#include "point.H"
#include "cellShapeControlMesh.H"
#include "triadField.H"
#include "scalarIOField.H"
#include "pointIOField.H"
#include "triadIOField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

template<class Triangulation, class Type>
Foam::tmp<Foam::Field<Type> > filterFarPoints
(
    const Triangulation& mesh,
    const Field<Type>& field
)
{
    tmp<Field<Type> > tNewField(new Field<Type>(field.size()));
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


template<class T>
autoPtr<mapDistribute> buildMap
(
    const T& mesh,
    labelListList& pointPoints
)
{
    pointPoints.setSize(mesh.vertexCount());

    globalIndex globalIndexing(mesh.vertexCount());

    for
    (
        typename T::Finite_vertices_iterator vit = mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        std::list<typename T::Vertex_handle> adjVerts;
        mesh.finite_adjacent_vertices(vit, std::back_inserter(adjVerts));

        DynamicList<label> indices(adjVerts.size());

        for
        (
            typename std::list<typename T::Vertex_handle>::const_iterator
                adjVertI = adjVerts.begin();
            adjVertI != adjVerts.end();
            ++adjVertI
        )
        {
            typename T::Vertex_handle vh = *adjVertI;

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

    List<Map<label> > compactMap;

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


template<class T>
Foam::tmp<Foam::triadField> buildAlignmentField(const T& mesh)
{
    tmp<triadField> tAlignments
    (
        new triadField(mesh.vertexCount(), triad::unset)
    );
    triadField& alignments = tAlignments();

    for
    (
        typename T::Finite_vertices_iterator vit = mesh.finite_vertices_begin();
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


template<class T>
Foam::tmp<Foam::pointField> buildPointField(const T& mesh)
{
    tmp<pointField> tPoints
    (
        new pointField(mesh.vertexCount(), point(GREAT, GREAT, GREAT))
    );
    pointField& points = tPoints();

    for
    (
        typename T::Finite_vertices_iterator vit = mesh.finite_vertices_begin();
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


void refine
(
    cellShapeControlMesh& mesh,
    const conformationSurfaces& geometryToConformTo,
    const label maxRefinementIterations,
    const scalar defaultCellSize
)
{
    for (label iter = 0; iter < maxRefinementIterations; ++iter)
    {
        DynamicList<point> ptsToInsert;

        for
        (
            CellSizeDelaunay::Finite_cells_iterator cit =
                mesh.finite_cells_begin();
            cit != mesh.finite_cells_end();
            ++cit
        )
        {
            const point newPoint =
                topoint
                (
                    CGAL::centroid
                    (
                        cit->vertex(0)->point(),
                        cit->vertex(1)->point(),
                        cit->vertex(2)->point(),
                        cit->vertex(3)->point()
                    )
                );

            if (geometryToConformTo.inside(newPoint))
            {
                ptsToInsert.append(newPoint);
            }
        }

        Info<< "    Adding " << returnReduce(ptsToInsert.size(), sumOp<label>())
            << endl;

        forAll(ptsToInsert, ptI)
        {
            mesh.insert
            (
                ptsToInsert[ptI],
                defaultCellSize,
                triad::unset,
                Vb::vtInternal
            );
        }
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    label maxRefinementIterations = 2;
    label maxSmoothingIterations = 200;
    scalar minResidual = 0;
    scalar defaultCellSize = 0.001;
    scalar nearFeatDistSqrCoeff = 1e-8;


    // Need to decouple vertex and cell type from this class?
    // Vertex must have:
    // + index
    // + procIndex
    // - type should be optional
    cellShapeControlMesh mesh(runTime);

    IOdictionary foamyHexMeshDict
    (
        IOobject
        (
            "foamyHexMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Random rndGen(64293*Pstream::myProcNo());

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "cvSearchableSurfaces",
            runTime.constant(),
            "triSurface",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        foamyHexMeshDict.subDict("geometry"),
        foamyHexMeshDict.lookupOrDefault("singleRegionName", true)
    );

    conformationSurfaces geometryToConformTo
    (
        runTime,
        rndGen,
        allGeometry,
        foamyHexMeshDict.subDict("surfaceConformation")
    );

    autoPtr<backgroundMeshDecomposition> bMesh;
    if (Pstream::parRun())
    {
        bMesh.set
        (
            new backgroundMeshDecomposition
            (
                runTime,
                rndGen,
                geometryToConformTo,
                foamyHexMeshDict.subDict("backgroundMeshDecomposition")
            )
        );
    }

    // Nice to have IO for the delaunay mesh
    // IO depend on vertex type.
    //
    // Define a delaunay mesh as:
    // + list of points of the triangulation
    // + optionally a list of cells

    Info<< nl << "Loop over surfaces" << endl;

    forAll(geometryToConformTo.surfaces(), sI)
    {
        const label surfI = geometryToConformTo.surfaces()[sI];

        const searchableSurface& surface =
            geometryToConformTo.geometry()[surfI];

        Info<< nl << "Inserting points from surface " << surface.name()
            << " (" << surface.type() << ")" << endl;

        const tmp<pointField> tpoints(surface.points());
        const pointField& points = tpoints();

        Info<< "    Number of points = " << points.size() << endl;

        forAll(points, pI)
        {
            // Is the point in the extendedFeatureEdgeMesh? If so get the
            // point normal, otherwise get the surface normal from
            // searchableSurface

            pointIndexHit info;
            label infoFeature;
            geometryToConformTo.findFeaturePointNearest
            (
                points[pI],
                nearFeatDistSqrCoeff,
                info,
                infoFeature
            );


            autoPtr<triad> pointAlignment;

            if (info.hit())
            {
                const extendedFeatureEdgeMesh& features =
                    geometryToConformTo.features()[infoFeature];

                vectorField norms = features.featurePointNormals(info.index());

                // Create a triad from these norms.
                pointAlignment.set(new triad());
                forAll(norms, nI)
                {
                    pointAlignment() += norms[nI];
                }

                pointAlignment().normalize();
                pointAlignment().orthogonalize();
            }
            else
            {
                geometryToConformTo.findEdgeNearest
                (
                    points[pI],
                    nearFeatDistSqrCoeff,
                    info,
                    infoFeature
                );

                if (info.hit())
                {
                    const extendedFeatureEdgeMesh& features =
                        geometryToConformTo.features()[infoFeature];

                    vectorField norms = features.edgeNormals(info.index());

                    // Create a triad from these norms.
                    pointAlignment.set(new triad());
                    forAll(norms, nI)
                    {
                        pointAlignment() += norms[nI];
                    }

                    pointAlignment().normalize();
                    pointAlignment().orthogonalize();
                }
                else
                {
                    pointField ptField(1, points[pI]);
                    scalarField distField(1, nearFeatDistSqrCoeff);
                    List<pointIndexHit> infoList(1, pointIndexHit());

                    surface.findNearest(ptField, distField, infoList);

                    vectorField normals(1);
                    surface.getNormal(infoList, normals);

                    pointAlignment.set(new triad(normals[0]));
                }
            }

            if (Pstream::parRun())
            {
                if (bMesh().positionOnThisProcessor(points[pI]))
                {
                    CellSizeDelaunay::Vertex_handle vh = mesh.insert
                    (
                        points[pI],
                        defaultCellSize,
                        pointAlignment(),
                        Vb::vtInternalNearBoundary
                    );
                }
            }
            else
            {
                CellSizeDelaunay::Vertex_handle vh = mesh.insert
                (
                    points[pI],
                    defaultCellSize,
                    pointAlignment(),
                    Vb::vtInternalNearBoundary
                );
            }
        }
    }


    // Refine the mesh
    refine
    (
        mesh,
        geometryToConformTo,
        maxRefinementIterations,
        defaultCellSize
    );


    if (Pstream::parRun())
    {
        mesh.distribute(bMesh);
    }


    labelListList pointPoints;
    autoPtr<mapDistribute> meshDistributor = buildMap(mesh, pointPoints);


    triadField alignments(buildAlignmentField(mesh));
    pointField points(buildPointField(mesh));

    mesh.printInfo(Info);


    // Setup the sizes and alignments on each point
    triadField fixedAlignments(mesh.vertexCount(), triad::unset);

    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (vit->nearBoundary())
        {
            fixedAlignments[vit->index()] = vit->alignment();
        }
    }

    Info<< nl << "Smoothing alignments" << endl;

    for (label iter = 0; iter < maxSmoothingIterations; iter++)
    {
        Info<< "Iteration " << iter;

        meshDistributor().distribute(points);
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

            const triad& oldTriad = alignments[pI];
            triad& newTriad = triadAv[pI];

            // Enforce the boundary conditions
            const triad& fixedAlignment = fixedAlignments[pI];

            forAll(pPoints, adjPointI)
            {
                const label adjPointIndex = pPoints[adjPointI];

                scalar dist = mag(points[pI] - points[adjPointIndex]);

//                dist = max(dist, SMALL);

                triad tmpTriad = alignments[adjPointIndex];

                for (direction dir = 0; dir < 3; dir++)
                {
                    if (tmpTriad.set(dir))
                    {
                        tmpTriad[dir] *= (1.0/dist);
                    }
                }

                newTriad += tmpTriad;
            }

            newTriad.normalize();
            newTriad.orthogonalize();
//            newTriad = newTriad.sortxyz();

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
                        newTriad.align(fixedAlignment[dirI]);
                    }
                }
            }
            else if (nFixed == 2)
            {
                forAll(fixedAlignment, dirI)
                {
                    if (fixedAlignment.set(dirI))
                    {
                        newTriad[dirI] = fixedAlignment[dirI];
                    }
                    else
                    {
                        newTriad[dirI] = triad::unset[dirI];
                    }
                }

                newTriad.orthogonalize();
            }
            else if (nFixed == 3)
            {
                forAll(fixedAlignment, dirI)
                {
                    if (fixedAlignment.set(dirI))
                    {
                        newTriad[dirI] = fixedAlignment[dirI];
                    }
                }
            }

            for (direction dir = 0; dir < 3; ++dir)
            {
                if
                (
                    newTriad.set(dir)
                 && oldTriad.set(dir)
                 //&& !fixedAlignment.set(dir)
                )
                {
                    scalar dotProd = (oldTriad[dir] & newTriad[dir]);
                    scalar diff = mag(dotProd) - 1.0;

                    residual += mag(diff);
                }
            }
        }

        forAll(alignments, pI)
        {
            alignments[pI] = triadAv[pI].sortxyz();
        }

        reduce(residual, sumOp<scalar>());

        Info<< ", Residual = " << residual << endl;

        if (residual <= minResidual)
        {
            break;
        }
    }


    // Write alignments to a .obj file
    OFstream str(runTime.path()/"alignments.obj");

    forAll(alignments, pI)
    {
        const triad& tri = alignments[pI];

        if (tri.set())
        {
            forAll(tri, dirI)
            {
                meshTools::writeOBJ(str, points[pI], tri[dirI] + points[pI]);
            }
        }
    }


    // Remove the far points
    pointIOField pointsIO
    (
        IOobject
        (
            "points",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        filterFarPoints(mesh, points)
    );

    scalarField sizes(points.size(), defaultCellSize);
    scalarIOField sizesIO
    (
        IOobject
        (
            "sizes",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        filterFarPoints(mesh, sizes)
    );

    triadIOField alignmentsIO
    (
        IOobject
        (
            "alignments",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        filterFarPoints(mesh, alignments)
    );

    pointsIO.write();
    sizesIO.write();
    alignmentsIO.write();

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
