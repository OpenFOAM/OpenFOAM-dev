/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "cutPolyIsoSurface.H"
#include "cellEdgeAddressing.H"
#include "cutPolyValue.H"
#include "EdgeMap.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cutPolyIsoSurface, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutPolyIsoSurface::cutPolyIsoSurface
(
    const polyMesh& mesh,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const labelList& zoneIDs
)
:
    points_(),
    pointEdges_(),
    pointEdgeLambdas_(),
    faces_(),
    faceCells_()
{
    cpuTime cpuTime;

    // Cut the faces
    List<List<labelPair>> faceCuts(mesh.faces().size());
    forAll(mesh.faces(), facei)
    {
        faceCuts[facei] =
            cutPoly::faceCuts
            (
                mesh.faces()[facei],
                pAlphas,
                isoAlpha
            );
    }

    // Request the cell-edge addressing engine
    const cellEdgeAddressingList& cAddrs = cellEdgeAddressingList::New(mesh);

    // Cut the cells
    List<labelListList> cellCuts(mesh.cells().size());
    label nCutCells = 0;
    auto cutCell = [&](const label celli)
    {
        cellCuts[celli] =
            cutPoly::cellCuts
            (
                mesh.cells()[celli],
                cAddrs[celli],
                mesh.faces(),
                faceCuts,
                pAlphas,
                isoAlpha
            );

        nCutCells += !cellCuts[celli].empty();
    };
    if (!isNull<labelList>(zoneIDs))
    {
        forAll(zoneIDs, i)
        {
            forAll(mesh.cellZones()[zoneIDs[i]], zoneCelli)
            {
                cutCell(mesh.cellZones()[zoneIDs[i]][zoneCelli]);
            }
        }
    }
    else
    {
        forAll(mesh.cells(), celli)
        {
            cutCell(celli);
        }
    }

    // Generate the cell cut polygons
    const label nAllocate = nCutCells + nCutCells/10;
    EdgeMap<label> meshEdgePoint(nAllocate*4);
    DynamicList<point> pointsDyn(nAllocate);
    DynamicList<edge> pointEdgesDyn(nAllocate);
    DynamicList<scalar> pointEdgeLambdasDyn(nAllocate);
    DynamicList<face> facesDyn(nAllocate);
    DynamicList<label> faceCellsDyn(nAllocate);
    forAll(mesh.cells(), celli)
    {
        forAll(cellCuts[celli], cellCuti)
        {
            if (cellCuts[celli][cellCuti].size() < 3) continue;

            facesDyn.append(face(cellCuts[celli][cellCuti].size()));

            forAll(cellCuts[celli][cellCuti], i)
            {
                const label cei = cellCuts[celli][cellCuti][i];
                const label cfi = cAddrs[celli].ceiToCfiAndFei()[cei][0][0];
                const label fei = cAddrs[celli].ceiToCfiAndFei()[cei][0][1];

                const edge e =
                    mesh.faces()[mesh.cells()[celli][cfi]].faceEdge(fei);

                EdgeMap<label>::iterator iter = meshEdgePoint.find(e);

                const label surfacePointi =
                    iter != meshEdgePoint.end() ? *iter : pointsDyn.size();

                if (surfacePointi == pointsDyn.size())
                {
                    meshEdgePoint.insert(e, surfacePointi);

                    const scalar lambda =
                        cutPoly::edgeCutLambda(e, pAlphas, isoAlpha);

                    pointsDyn.append
                    (
                        cutPoly::edgeCutValue(e, lambda, mesh.points())
                    );
                    pointEdgesDyn.append(e);
                    pointEdgeLambdasDyn.append(lambda);
                }

                facesDyn.last()[i] = surfacePointi;
            }

            faceCellsDyn.append(celli);
        }
    }

    // Transfer to non-dynamic storage
    points_.transfer(pointsDyn);
    pointEdges_.transfer(pointEdgesDyn);
    pointEdgeLambdas_.transfer(pointEdgeLambdasDyn);
    faces_.transfer(facesDyn);
    faceCells_.transfer(faceCellsDyn);

    if (debug)
    {
        const label nLabels =
            sum(ListListOps::subSizes(faces(), accessOp<face>()));

        Pout<< typeName << " : constructed surface of size "
            << points().size()*sizeof(point) + nLabels*sizeof(label)
            << " in " << cpuTime.cpuTimeIncrement() << "s" << nl << endl;
    }
}


Foam::cutPolyIsoSurface::cutPolyIsoSurface
(
    const PtrList<cutPolyIsoSurface>& isos
)
:
    points_(),
    pointEdges_(),
    pointEdgeLambdas_(),
    faces_(),
    faceCells_()
{
    label nPoints = 0, nFaces = 0;
    forAll(isos, i)
    {
        nPoints += isos[i].points_.size();
        nFaces += isos[i].faces_.size();
    }

    points_.resize(nPoints);
    pointEdges_.resize(nPoints);
    pointEdgeLambdas_.resize(nPoints);
    faces_.resize(nFaces);
    faceCells_.resize(nFaces);

    label pointi0 = 0, facei0 = 0;
    forAll(isos, i)
    {
        const label nPs = isos[i].points_.size();
        const label nFs = isos[i].faces_.size();

        SubList<point>(points_, nPs, pointi0) = isos[i].points_;
        SubList<edge>(pointEdges_, nPs, pointi0) = isos[i].pointEdges_;
        SubList<scalar>(pointEdgeLambdas_, nPs, pointi0) =
            isos[i].pointEdgeLambdas_;
        SubList<face>(faces_, nFs, facei0) = isos[i].faces_;
        SubList<label>(faceCells_, nFs, facei0) = isos[i].faceCells_;

        forAll(isos[i].faces_, fi)
        {
            forAll(faces_[facei0 + fi], fpi)
            {
                faces_[facei0 + fi][fpi] += pointi0;
            }
        }

        pointi0 += nPs;
        facei0 += nFs;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cutPolyIsoSurface::~cutPolyIsoSurface()
{}


// ************************************************************************* //
