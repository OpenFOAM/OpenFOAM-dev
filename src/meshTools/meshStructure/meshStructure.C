/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    const label cellI
) const
{
    const cell& cFaces = mesh.cells()[cellI];

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
                label pointI = f[fp];
                if (pointLayer_[pointI] == layerI)
                {
                    nLayer++;
                }
                else if (pointLayer_[pointI] == layerI+1)
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
        forAll(pp, patchFaceI)
        {
            patchFaces[patchFaceI] = pp.addressing()[patchFaceI];
            patchData[patchFaceI] = topoDistanceData(patchFaceI, 0);
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
        forAll(cellToPatchFaceAddressing_, cellI)
        {
            cellToPatchFaceAddressing_[cellI] = cellData[cellI].data();
            cellLayer_[cellI] = cellData[cellI].distance();
        }



        // Determine faces from face-cell-walk
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        faceToPatchFaceAddressing_.setSize(mesh.nFaces());
        faceToPatchEdgeAddressing_.setSize(mesh.nFaces());
        faceToPatchEdgeAddressing_ = labelMin;
        faceLayer_.setSize(mesh.nFaces());

        forAll(faceToPatchFaceAddressing_, faceI)
        {
            label own = mesh.faceOwner()[faceI];
            label patchFaceI = faceData[faceI].data();
            label patchDist = faceData[faceI].distance();

            if (mesh.isInternalFace(faceI))
            {
                label nei = mesh.faceNeighbour()[faceI];

                if (cellData[own].distance() == cellData[nei].distance())
                {
                    // side face
                    faceToPatchFaceAddressing_[faceI] = 0;
                    faceLayer_[faceI] = cellData[own].distance();
                }
                else if (cellData[own].distance() < cellData[nei].distance())
                {
                    // unturned face
                    faceToPatchFaceAddressing_[faceI] = patchFaceI+1;
                    faceToPatchEdgeAddressing_[faceI] = -1;
                    faceLayer_[faceI] = patchDist;
                }
                else
                {
                    // turned face
                    faceToPatchFaceAddressing_[faceI] = -(patchFaceI+1);
                    faceToPatchEdgeAddressing_[faceI] = -1;
                    faceLayer_[faceI] = patchDist;
                }
            }
            else if (patchDist == cellData[own].distance())
            {
                // starting face
                faceToPatchFaceAddressing_[faceI] = -(patchFaceI+1);
                faceToPatchEdgeAddressing_[faceI] = -1;
                faceLayer_[faceI] = patchDist;
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
        forAll(pp.meshPoints(), patchPointI)
        {
            patchPoints[patchPointI] = pp.meshPoints()[patchPointI];
            patchData[patchPointI] = pointTopoDistanceData(patchPointI, 0);
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

        forAll(pointData, pointI)
        {
            pointToPatchPointAddressing_[pointI] = pointData[pointI].data();
            pointLayer_[pointI] = pointData[pointI].distance();
        }


        // Derive from originating patch points what the patch edges were.
        EdgeMap<label> pointsToEdge(pp.nEdges());
        forAll(pp.edges(), edgeI)
        {
            pointsToEdge.insert(pp.edges()[edgeI], edgeI);
        }

        // Look up on faces
        forAll(faceToPatchEdgeAddressing_, faceI)
        {
            if (faceToPatchEdgeAddressing_[faceI] == labelMin)
            {
                // Face not yet done. Check if all points on same level
                // or if not see what edge it originates from

                const face& f = mesh.faces()[faceI];

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
                    //Pout<< "Horizontal boundary face " << faceI
                    //    << " at:" << mesh.faceCentres()[faceI]
                    //    << " data:" << faceData[faceI]
                    //    << " pointDatas:"
                    //    << UIndirectList<pointTopoDistanceData>(pointData, f)
                    //    << endl;

                    label patchFaceI = faceData[faceI].data();
                    label patchDist = faceData[faceI].distance();

                    faceToPatchEdgeAddressing_[faceI] = -1;
                    faceToPatchFaceAddressing_[faceI] = patchFaceI+1;
                    faceLayer_[faceI] = patchDist;
                }
                else
                {
                    // Points of face on different levels

                    // See if there is any edge
                    forAll(f, fp)
                    {
                        label pointI = f[fp];
                        label nextPointI = f.nextLabel(fp);

                        EdgeMap<label>::const_iterator fnd = pointsToEdge.find
                        (
                            edge
                            (
                                pointData[pointI].data(),
                                pointData[nextPointI].data()
                            )
                        );
                        if (fnd != pointsToEdge.end())
                        {
                            faceToPatchEdgeAddressing_[faceI] = fnd();
                            faceToPatchFaceAddressing_[faceI] = 0;
                            label own = mesh.faceOwner()[faceI];
                            faceLayer_[faceI] = cellData[own].distance();

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

            forAll(lCells, lCellI)
            {
                label cellI = lCells[lCellI];

                structured_ = isStructuredCell
                (
                    mesh,
                    layerI,
                    cellI
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
