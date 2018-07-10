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

#include "meshStructure.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "pointTopoDistanceData.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(meshStructure, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshStructure::isStructuredCell
(
    const polyMesh& mesh,
    const label layerI,
    const label celli
) const
{
    const cell& cFaces = mesh.cells()[celli];

    // Count number of side faces
    label nSide = 0;
    forAll(cFaces, i)
    {
        if (faceToPatchEdgeAddressing_[cFaces[i]] != -1)
        {
            nSide++;
        }
    }

    if (nSide != cFaces.size()-2)
    {
        return false;
    }

    // Check that side faces have correct point layers
    forAll(cFaces, i)
    {
        if (faceToPatchEdgeAddressing_[cFaces[i]] != -1)
        {
            const face& f = mesh.faces()[cFaces[i]];

            label nLayer = 0;
            label nLayerPlus1 = 0;
            forAll(f, fp)
            {
                label pointi = f[fp];
                if (pointLayer_[pointi] == layerI)
                {
                    nLayer++;
                }
                else if (pointLayer_[pointi] == layerI+1)
                {
                    nLayerPlus1++;
                }
            }

            if (f.size() != 4 || (nLayer+nLayerPlus1 != 4))
            {
                return false;
            }
        }
    }

    return true;
}


void Foam::meshStructure::correct
(
    const polyMesh& mesh,
    const uindirectPrimitivePatch& pp
)
{
    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    {
        if (debug)
        {
            Info<< typeName << " : seeding "
                << returnReduce(pp.size(), sumOp<label>()) << " patch faces"
                << nl << endl;
        }


        // Start of changes
        labelList patchFaces(pp.size());
        List<topoDistanceData> patchData(pp.size());
        forAll(pp, patchFacei)
        {
            patchFaces[patchFacei] = pp.addressing()[patchFacei];
            patchData[patchFacei] = topoDistanceData(patchFacei, 0);
        }


        // Propagate information inwards
        FaceCellWave<topoDistanceData> distanceCalc
        (
            mesh,
            patchFaces,
            patchData,
            faceData,
            cellData,
            mesh.globalData().nTotalCells()+1
        );


        // Determine cells from face-cell-walk
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        cellToPatchFaceAddressing_.setSize(mesh.nCells());
        cellLayer_.setSize(mesh.nCells());
        forAll(cellToPatchFaceAddressing_, celli)
        {
            cellToPatchFaceAddressing_[celli] = cellData[celli].data();
            cellLayer_[celli] = cellData[celli].distance();
        }



        // Determine faces from face-cell-walk
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        faceToPatchFaceAddressing_.setSize(mesh.nFaces());
        faceToPatchEdgeAddressing_.setSize(mesh.nFaces());
        faceToPatchEdgeAddressing_ = labelMin;
        faceLayer_.setSize(mesh.nFaces());

        forAll(faceToPatchFaceAddressing_, facei)
        {
            label own = mesh.faceOwner()[facei];
            label patchFacei = faceData[facei].data();
            label patchDist = faceData[facei].distance();

            if (mesh.isInternalFace(facei))
            {
                label nei = mesh.faceNeighbour()[facei];

                if (cellData[own].distance() == cellData[nei].distance())
                {
                    // side face
                    faceToPatchFaceAddressing_[facei] = 0;
                    faceLayer_[facei] = cellData[own].distance();
                }
                else if (cellData[own].distance() < cellData[nei].distance())
                {
                    // unturned face
                    faceToPatchFaceAddressing_[facei] = patchFacei+1;
                    faceToPatchEdgeAddressing_[facei] = -1;
                    faceLayer_[facei] = patchDist;
                }
                else
                {
                    // turned face
                    faceToPatchFaceAddressing_[facei] = -(patchFacei+1);
                    faceToPatchEdgeAddressing_[facei] = -1;
                    faceLayer_[facei] = patchDist;
                }
            }
            else if (patchDist == cellData[own].distance())
            {
                // starting face
                faceToPatchFaceAddressing_[facei] = -(patchFacei+1);
                faceToPatchEdgeAddressing_[facei] = -1;
                faceLayer_[facei] = patchDist;
            }
            else
            {
                // unturned face or side face. Cannot be determined until
                // we determine the point layers. Problem is that both are
                // the same number of steps away from the initial seed face.
            }
        }
    }


    // Determine points from separate walk on point-edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        pointToPatchPointAddressing_.setSize(mesh.nPoints());
        pointLayer_.setSize(mesh.nPoints());

        if (debug)
        {
            Info<< typeName << " : seeding "
                << returnReduce(pp.nPoints(), sumOp<label>()) << " patch points"
                << nl << endl;
        }

        // Field on edges and points.
        List<pointTopoDistanceData> edgeData(mesh.nEdges());
        List<pointTopoDistanceData> pointData(mesh.nPoints());

        // Start of changes
        labelList patchPoints(pp.nPoints());
        List<pointTopoDistanceData> patchData(pp.nPoints());
        forAll(pp.meshPoints(), patchPointi)
        {
            patchPoints[patchPointi] = pp.meshPoints()[patchPointi];
            patchData[patchPointi] = pointTopoDistanceData(patchPointi, 0);
        }


        // Walk
        PointEdgeWave<pointTopoDistanceData> distanceCalc
        (
            mesh,
            patchPoints,
            patchData,

            pointData,
            edgeData,
            mesh.globalData().nTotalPoints()  // max iterations
        );

        forAll(pointData, pointi)
        {
            pointToPatchPointAddressing_[pointi] = pointData[pointi].data();
            pointLayer_[pointi] = pointData[pointi].distance();
        }


        // Derive from originating patch points what the patch edges were.
        EdgeMap<label> pointsToEdge(pp.nEdges());
        forAll(pp.edges(), edgeI)
        {
            pointsToEdge.insert(pp.edges()[edgeI], edgeI);
        }

        // Look up on faces
        forAll(faceToPatchEdgeAddressing_, facei)
        {
            if (faceToPatchEdgeAddressing_[facei] == labelMin)
            {
                // Face not yet done. Check if all points on same level
                // or if not see what edge it originates from

                const face& f = mesh.faces()[facei];

                label levelI = pointLayer_[f[0]];
                for (label fp = 1; fp < f.size(); fp++)
                {
                    if (pointLayer_[f[fp]] != levelI)
                    {
                        levelI = -1;
                        break;
                    }
                }

                if (levelI != -1)
                {
                    // All same level
                    // Pout<< "Horizontal boundary face " << facei
                    //    << " at:" << mesh.faceCentres()[facei]
                    //    << " data:" << faceData[facei]
                    //    << " pointDatas:"
                    //    << UIndirectList<pointTopoDistanceData>(pointData, f)
                    //    << endl;

                    label patchFacei = faceData[facei].data();
                    label patchDist = faceData[facei].distance();

                    faceToPatchEdgeAddressing_[facei] = -1;
                    faceToPatchFaceAddressing_[facei] = patchFacei+1;
                    faceLayer_[facei] = patchDist;
                }
                else
                {
                    // Points of face on different levels

                    // See if there is any edge
                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        label nextPointi = f.nextLabel(fp);

                        EdgeMap<label>::const_iterator fnd = pointsToEdge.find
                        (
                            edge
                            (
                                pointData[pointi].data(),
                                pointData[nextPointi].data()
                            )
                        );
                        if (fnd != pointsToEdge.end())
                        {
                            faceToPatchEdgeAddressing_[facei] = fnd();
                            faceToPatchFaceAddressing_[facei] = 0;
                            label own = mesh.faceOwner()[facei];
                            faceLayer_[facei] = cellData[own].distance();

                            // Note: could test whether the other edges on the
                            // face are consistent
                            break;
                        }
                    }
                }
            }
        }
    }



    // Use maps to find out mesh structure.
    {
        label nLayers = gMax(cellLayer_)+1;
        labelListList layerToCells(invertOneToMany(nLayers, cellLayer_));

        structured_ = true;
        forAll(layerToCells, layerI)
        {
            const labelList& lCells = layerToCells[layerI];

            forAll(lCells, lCelli)
            {
                label celli = lCells[lCelli];

                structured_ = isStructuredCell
                (
                    mesh,
                    layerI,
                    celli
                );

                if (!structured_)
                {
                    break;
                }
            }

            if (!structured_)
            {
                break;
            }
        }

        reduce(structured_, andOp<bool>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshStructure::meshStructure
(
    const polyMesh& mesh,
    const uindirectPrimitivePatch& pp
)
{
    correct(mesh, pp);
}


// ************************************************************************* //
