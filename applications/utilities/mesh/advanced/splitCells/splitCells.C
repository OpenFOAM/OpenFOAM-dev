/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    splitCells

Description
    Utility to split cells with flat faces.

    Uses a geometric cut with a plane dividing the edge angle into two so
    might produce funny cells. For hexes it will use by default a cut from
    edge onto opposite edge (i.e. purely topological).

    Options:
    - split cells from cellSet only
    - use geometric cut for hexes as well

    The angle is the angle between two faces sharing an edge as seen from
    inside each cell. So a cube will have all angles 90. If you want
    to split cells with cell centre outside use e.g. angle 200

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "cellCuts.H"
#include "cellSet.H"
#include "cellModeller.H"
#include "meshCutter.H"
#include "unitConversion.H"
#include "geomCellLooper.H"
#include "plane.H"
#include "edgeVertex.H"
#include "meshTools.H"
#include "ListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


labelList pack(const boolList& lst)
{
    labelList packedLst(lst.size());
    label packedI = 0;

    forAll(lst, i)
    {
        if (lst[i])
        {
            packedLst[packedI++] = i;
        }
    }
    packedLst.setSize(packedI);

    return packedLst;
}


scalarField pack(const boolList& lst, const scalarField& elems)
{
    scalarField packedElems(lst.size());
    label packedI = 0;

    forAll(lst, i)
    {
        if (lst[i])
        {
            packedElems[packedI++] = elems[i];
        }
    }
    packedElems.setSize(packedI);

    return packedElems;
}


// Given sin and cos of max angle between normals calculate whether f0 and f1
// on celli make larger angle. Uses sinAngle only for quadrant detection.
bool largerAngle
(
    const primitiveMesh& mesh,
    const scalar cosAngle,
    const scalar sinAngle,

    const label celli,
    const label f0,             // face label
    const label f1,

    const vector& n0,           // normal at f0
    const vector& n1
)
{
    const labelList& own = mesh.faceOwner();

    bool sameFaceOrder = !((own[f0] == celli) ^ (own[f1] == celli));

    // Get cos between faceArea vectors. Correct so flat angle (180 degrees)
    // gives -1.
    scalar normalCosAngle = n0 & n1;

    if (sameFaceOrder)
    {
        normalCosAngle = -normalCosAngle;
    }


    // Get cos between faceCentre and normal vector to determine in
    // which quadrant angle is. (Is correct for unwarped faces only!)
    // Correct for non-outwards pointing normal.
    vector c1c0(mesh.faceCentres()[f1] - mesh.faceCentres()[f0]);
    c1c0 /= mag(c1c0) + vSmall;

    scalar fcCosAngle = n0 & c1c0;

    if (own[f0] != celli)
    {
        fcCosAngle = -fcCosAngle;
    }

    if (sinAngle < 0.0)
    {
        // Looking for concave angles (quadrant 3 or 4)
        if (fcCosAngle <= 0)
        {
            // Angle is convex so smaller.
            return false;
        }
        else
        {
            if (normalCosAngle < cosAngle)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
    else
    {
        // Looking for convex angles (quadrant 1 or 2)
        if (fcCosAngle > 0)
        {
            // Concave angle
            return true;
        }
        else
        {
            // Convex. Check cos of normal vectors.
            if (normalCosAngle > cosAngle)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
}


// Split hex (and hex only) along edgeI creating two prisms
bool splitHex
(
    const polyMesh& mesh,
    const label celli,
    const label edgeI,

    DynamicList<label>& cutCells,
    DynamicList<labelList>& cellLoops,
    DynamicList<scalarField>& cellEdgeWeights
)
{
    // cut handling functions
    edgeVertex ev(mesh);

    const edgeList& edges = mesh.edges();
    const faceList& faces = mesh.faces();

    const edge& e = edges[edgeI];

    // Get faces on the side, i.e. faces not using edge but still using one of
    // the edge endpoints.

    label leftI = -1;
    label rightI = -1;
    label leftFp = -1;
    label rightFp = -1;

    const cell& cFaces = mesh.cells()[celli];

    forAll(cFaces, i)
    {
        label facei = cFaces[i];

        const face& f = faces[facei];

        label fp0 = findIndex(f, e[0]);
        label fp1 = findIndex(f, e[1]);

        if (fp0 == -1)
        {
            if (fp1 != -1)
            {
                // Face uses e[1] but not e[0]
                rightI = facei;
                rightFp = fp1;

                if (leftI != -1)
                {
                    // Have both faces so exit
                    break;
                }
            }
        }
        else
        {
            if (fp1 != -1)
            {
                // Face uses both e[1] and e[0]
            }
            else
            {
                leftI = facei;
                leftFp = fp0;

                if (rightI != -1)
                {
                    break;
                }
            }
        }
    }

    if (leftI == -1 || rightI == -1)
    {
        FatalErrorInFunction
            << " rightI:" << rightI << abort(FatalError);
    }

    // Walk two vertices further on faces.

    const face& leftF = faces[leftI];

    label leftV = leftF[(leftFp + 2) % leftF.size()];

    const face& rightF = faces[rightI];

    label rightV = rightF[(rightFp + 2) % rightF.size()];

    labelList loop(4);
    loop[0] = ev.vertToEVert(e[0]);
    loop[1] = ev.vertToEVert(leftV);
    loop[2] = ev.vertToEVert(rightV);
    loop[3] = ev.vertToEVert(e[1]);

    scalarField loopWeights(4);
    loopWeights[0] = -great;
    loopWeights[1] = -great;
    loopWeights[2] = -great;
    loopWeights[3] = -great;

    cutCells.append(celli);
    cellLoops.append(loop);
    cellEdgeWeights.append(loopWeights);

    return true;
}


// Split celli along edgeI with a plane along halfNorm direction.
bool splitCell
(
    const polyMesh& mesh,
    const geomCellLooper& cellCutter,

    const label celli,
    const label edgeI,
    const vector& halfNorm,

    const boolList& vertIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,

    DynamicList<label>& cutCells,
    DynamicList<labelList>& cellLoops,
    DynamicList<scalarField>& cellEdgeWeights
)
{
    const edge& e = mesh.edges()[edgeI];

    vector eVec = e.vec(mesh.points());
    eVec /= mag(eVec);

    vector planeN = eVec ^ halfNorm;

    // Slightly tilt plane to make it not cut edges exactly
    // halfway on fully regular meshes (since we want cuts
    // to be snapped to vertices)
    planeN += 0.01*halfNorm;

    planeN /= mag(planeN);

    // Define plane through edge
    plane cutPlane(mesh.points()[e.start()], planeN);

    labelList loop;
    scalarField loopWeights;

    if
    (
        cellCutter.cut
        (
            cutPlane,
            celli,
            vertIsCut,
            edgeIsCut,
            edgeWeight,
            loop,
            loopWeights
        )
    )
    {
        // Did manage to cut cell. Copy into overall list.
        cutCells.append(celli);
        cellLoops.append(loop);
        cellEdgeWeights.append(loopWeights);

        return true;
    }
    else
    {
        return false;
    }
}


// Collects cuts for all cells in cellSet
void collectCuts
(
    const polyMesh& mesh,
    const geomCellLooper& cellCutter,
    const bool geometry,
    const scalar minCos,
    const scalar minSin,
    const cellSet& cellsToCut,

    DynamicList<label>& cutCells,
    DynamicList<labelList>& cellLoops,
    DynamicList<scalarField>& cellEdgeWeights
)
{
    // Get data from mesh
    const cellShapeList& cellShapes = mesh.cellShapes();
    const labelList& own = mesh.faceOwner();
    const labelListList& cellEdges = mesh.cellEdges();
    const vectorField& faceAreas = mesh.faceAreas();

    // Hex shape
    const cellModel& hex = *(cellModeller::lookup("hex"));

    // cut handling functions
    edgeVertex ev(mesh);


    // Cut information per mesh entity
    boolList vertIsCut(mesh.nPoints(), false);
    boolList edgeIsCut(mesh.nEdges(), false);
    scalarField edgeWeight(mesh.nEdges(), -great);

    forAllConstIter(cellSet, cellsToCut, iter)
    {
        const label celli = iter.key();
        const labelList& cEdges = cellEdges[celli];

        forAll(cEdges, i)
        {
            label edgeI = cEdges[i];

            label f0, f1;
            meshTools::getEdgeFaces(mesh, celli, edgeI, f0, f1);

            vector n0 = faceAreas[f0];
            n0 /= mag(n0);

            vector n1 = faceAreas[f1];
            n1 /= mag(n1);

            if
            (
                largerAngle
                (
                    mesh,
                    minCos,
                    minSin,

                    celli,
                    f0,
                    f1,
                    n0,
                    n1
                )
            )
            {
                bool splitOk = false;

                if (!geometry && cellShapes[celli].model() == hex)
                {
                    splitOk =
                        splitHex
                        (
                            mesh,
                            celli,
                            edgeI,

                            cutCells,
                            cellLoops,
                            cellEdgeWeights
                        );
                }
                else
                {
                    vector halfNorm;

                    if ((own[f0] == celli) ^ (own[f1] == celli))
                    {
                        // Opposite owner orientation
                        halfNorm = 0.5*(n0 - n1);
                    }
                    else
                    {
                        // Faces have same owner or same neighbour so
                        // normals point in same direction
                        halfNorm = 0.5*(n0 + n1);
                    }

                    splitOk =
                        splitCell
                        (
                            mesh,
                            cellCutter,
                            celli,
                            edgeI,
                            halfNorm,

                            vertIsCut,
                            edgeIsCut,
                            edgeWeight,

                            cutCells,
                            cellLoops,
                            cellEdgeWeights
                        );
                }


                if (splitOk)
                {
                    // Update cell/edge/vertex wise info.
                    label index = cellLoops.size() - 1;
                    const labelList& loop = cellLoops[index];
                    const scalarField& loopWeights = cellEdgeWeights[index];

                    forAll(loop, i)
                    {
                        label cut = loop[i];

                        if (ev.isEdge(cut))
                        {
                            edgeIsCut[ev.getEdge(cut)] = true;
                            edgeWeight[ev.getEdge(cut)] = loopWeights[i];
                        }
                        else
                        {
                            vertIsCut[ev.getVertex(cut)] = true;
                        }
                    }

                    // Stop checking edges for this cell.
                    break;
                }
            }
        }
    }

    cutCells.shrink();
    cellLoops.shrink();
    cellEdgeWeights.shrink();
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "split cells with flat faces"
    );
    #include "addOverwriteOption.H"
    argList::noParallel();
    argList::validArgs.append("edgeAngle [0..360]");

    argList::addOption
    (
        "set",
        "name",
        "split cells from specified cellSet only"
    );
    argList::addBoolOption
    (
        "geometry",
        "use geometric cut for hexes as well"
    );
    argList::addOption
    (
        "tol",
        "scalar", "edge snap tolerance (default 0.2)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const scalar featureAngle = args.argRead<scalar>(1);
    const scalar minCos = Foam::cos(degToRad(featureAngle));
    const scalar minSin = Foam::sin(degToRad(featureAngle));

    const bool readSet   = args.optionFound("set");
    const bool geometry  = args.optionFound("geometry");
    const bool overwrite = args.optionFound("overwrite");

    const scalar edgeTol = args.optionLookupOrDefault("tol", 0.2);

    Info<< "Trying to split cells with internal angles > feature angle\n" << nl
        << "featureAngle      : " << featureAngle << nl
        << "edge snapping tol : " << edgeTol << nl;
    if (readSet)
    {
        Info<< "candidate cells   : cellSet " << args["set"] << nl;
    }
    else
    {
        Info<< "candidate cells   : all cells" << nl;
    }
    if (geometry)
    {
        Info<< "hex cuts          : geometric; using edge tolerance" << nl;
    }
    else
    {
        Info<< "hex cuts          : topological; cut to opposite edge" << nl;
    }
    Info<< endl;


    // Cell circumference cutter
    geomCellLooper cellCutter(mesh);
    // Snap all edge cuts close to endpoints to vertices.
    geomCellLooper::setSnapTol(edgeTol);

    // Candidate cells to cut
    cellSet cellsToCut(mesh, "cellsToCut", mesh.nCells()/100);

    if (readSet)
    {
        // Read cells to cut from cellSet
        cellSet cells(mesh, args["set"]);

        cellsToCut = cells;
    }

    while (true)
    {
        if (!readSet)
        {
            // Try all cells for cutting
            for (label celli = 0; celli < mesh.nCells(); celli++)
            {
                cellsToCut.insert(celli);
            }
        }


        // Cut information per cut cell
        DynamicList<label> cutCells(mesh.nCells()/10 + 10);
        DynamicList<labelList> cellLoops(mesh.nCells()/10 + 10);
        DynamicList<scalarField> cellEdgeWeights(mesh.nCells()/10 + 10);

        collectCuts
        (
            mesh,
            cellCutter,
            geometry,
            minCos,
            minSin,
            cellsToCut,

            cutCells,
            cellLoops,
            cellEdgeWeights
        );

        cellSet cutSet(mesh, "cutSet", cutCells.size());
        forAll(cutCells, i)
        {
            cutSet.insert(cutCells[i]);
        }

        // Gets cuts across cells from cuts through edges.
        Info<< "Writing " << cutSet.size() << " cells to cut to cellSet "
            << cutSet.instance()/cutSet.local()/cutSet.name()
            << endl << endl;
        cutSet.write();

        // Analyze cuts for clashes.
        cellCuts cuts
        (
            mesh,
            cutCells,       // cells candidate for cutting
            cellLoops,
            cellEdgeWeights
        );

        Info<< "Actually cut cells:" << cuts.nLoops() << nl << endl;

        if (cuts.nLoops() == 0)
        {
            break;
        }

        // Remove cut cells from cellsToCut  (Note:only relevant if -readSet)
        forAll(cuts.cellLoops(), celli)
        {
            if (cuts.cellLoops()[celli].size())
            {
                // Info<< "Removing cut cell " << celli << " from wishlist"
                //    << endl;
                cellsToCut.erase(celli);
            }
        }

        // At least some cells are cut.
        polyTopoChange meshMod(mesh);

        // Cutting engine
        meshCutter cutter(mesh);

        // Insert mesh refinement into polyTopoChange.
        cutter.setRefinement(cuts, meshMod);

        // Do all changes
        Info<< "Morphing ..." << endl;

        if (!overwrite)
        {
            runTime++;
        }

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Update stored labels on meshCutter
        cutter.updateMesh(morphMap());

        // Update cellSet
        cellsToCut.updateMesh(morphMap());

        Info<< "Remaining:" << cellsToCut.size() << endl;

        // Write resulting mesh
        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }

        Info<< "Writing refined morphMesh to time " << runTime.timeName()
            << endl;

        mesh.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
