/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    surfaceFeatureExtract

Description
    Extracts and writes surface features to file. All but the basic feature
    extraction is WIP.

    Curvature calculation is an implementation of the algorithm from:

      "Estimating Curvatures and their Derivatives on Triangle Meshes"
      by S. Rusinkiewicz

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface.H"
#include "surfaceFeatures.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "triSurfaceFields.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "unitConversion.H"
#include "plane.H"
#include "tensor2D.H"
#include "symmTensor2D.H"
#include "point.H"
#include "triadField.H"
#include "transform.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar calcVertexNormalWeight
(
    const triFace& f,
    const label pI,
    const vector& fN,
    const pointField& points
)
{
    label index = findIndex(f, pI);

    if (index == -1)
    {
        FatalErrorIn("calcVertexNormals()")
            << "Point not in face" << abort(FatalError);
    }

    const vector e1 = points[f[index]] - points[f[f.fcIndex(index)]];
    const vector e2 = points[f[index]] - points[f[f.rcIndex(index)]];

    return mag(fN)/(magSqr(e1)*magSqr(e2) + VSMALL);
}


point randomPointInPlane(const plane& p)
{
    // Perturb base point
    const point& refPt = p.refPoint();

    // ax + by + cz + d = 0
    const FixedList<scalar, 4>& planeCoeffs = p.planeCoeffs();

    const scalar perturbX = refPt.x() + 1e-3;
    const scalar perturbY = refPt.y() + 1e-3;
    const scalar perturbZ = refPt.z() + 1e-3;

    if (mag(planeCoeffs[2]) < SMALL)
    {
        if (mag(planeCoeffs[1]) < SMALL)
        {
            const scalar x =
                -1.0
                *(
                     planeCoeffs[3]
                   + planeCoeffs[1]*perturbY
                   + planeCoeffs[2]*perturbZ
                 )/planeCoeffs[0];

            return point
            (
                x,
                perturbY,
                perturbZ
            );
        }

        const scalar y =
            -1.0
            *(
                 planeCoeffs[3]
               + planeCoeffs[0]*perturbX
               + planeCoeffs[2]*perturbZ
             )/planeCoeffs[1];

        return point
        (
            perturbX,
            y,
            perturbZ
        );
    }
    else
    {
        const scalar z =
            -1.0
            *(
                 planeCoeffs[3]
               + planeCoeffs[0]*perturbX
               + planeCoeffs[1]*perturbY
             )/planeCoeffs[2];

        return point
        (
            perturbX,
            perturbY,
            z
        );
    }
}


triadField calcVertexCoordSys
(
    const triSurface& surf,
    const vectorField& pointNormals
)
{
    const pointField& points = surf.points();
    const Map<label>& meshPointMap = surf.meshPointMap();

    triadField pointCoordSys(points.size());

    forAll(points, pI)
    {
        const point& pt = points[pI];
        const vector& normal = pointNormals[meshPointMap[pI]];

        if (mag(normal) < SMALL)
        {
            pointCoordSys[meshPointMap[pI]] = triad::unset;
            continue;
        }

        plane p(pt, normal);

        // Pick random point in plane
        vector dir1 = pt - randomPointInPlane(p);
        dir1 /= mag(dir1);

        vector dir2 = dir1 ^ normal;
        dir2 /= mag(dir2);

        pointCoordSys[meshPointMap[pI]] = triad(dir1, dir2, normal);
    }

    return pointCoordSys;
}


vectorField calcVertexNormals(const triSurface& surf)
{
    // Weighted average of normals of faces attached to the vertex
    // Weight = fA / (mag(e0)^2 * mag(e1)^2);

    Info<< "Calculating vertex normals" << endl;

    vectorField pointNormals(surf.nPoints(), vector::zero);

    const pointField& points = surf.points();
    const labelListList& pointFaces = surf.pointFaces();
    const labelList& meshPoints = surf.meshPoints();

    forAll(pointFaces, pI)
    {
        const labelList& pFaces = pointFaces[pI];

        forAll(pFaces, fI)
        {
            const label faceI = pFaces[fI];
            const triFace& f = surf[faceI];

            vector fN = f.normal(points);

            scalar weight = calcVertexNormalWeight
            (
                f,
                meshPoints[pI],
                fN,
                points
            );

            pointNormals[pI] += weight*fN;
        }

        pointNormals[pI] /= mag(pointNormals[pI]) + VSMALL;
    }

    return pointNormals;
}


triSurfacePointScalarField calcCurvature
(
    const word& name,
    const Time& runTime,
    const triSurface& surf,
    const vectorField& pointNormals,
    const triadField& pointCoordSys
)
{
    Info<< "Calculating face curvature" << endl;

    const pointField& points = surf.points();
    const labelList& meshPoints = surf.meshPoints();
    const Map<label>& meshPointMap = surf.meshPointMap();

    triSurfacePointScalarField curvaturePointField
    (
        IOobject
        (
            name + ".curvature",
            runTime.constant(),
            "triSurface",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        scalarField(points.size(), 0.0)
    );

    List<symmTensor2D> pointFundamentalTensors
    (
        points.size(),
        symmTensor2D::zero
    );

    scalarList accumulatedWeights(points.size(), 0.0);

    forAll(surf, fI)
    {
        const triFace& f = surf[fI];
        const edgeList fEdges = f.edges();

        // Calculate the edge vectors and the normal differences
        vectorField edgeVectors(f.size(), vector::zero);
        vectorField normalDifferences(f.size(), vector::zero);

        forAll(fEdges, feI)
        {
            const edge& e = fEdges[feI];

            edgeVectors[feI] = e.vec(points);
            normalDifferences[feI] =
               pointNormals[meshPointMap[e[0]]]
             - pointNormals[meshPointMap[e[1]]];
        }

        // Set up a local coordinate system for the face
        const vector& e0 = edgeVectors[0];
        const vector eN = f.normal(points);
        const vector e1 = (e0 ^ eN);

        if (magSqr(eN) < ROOTVSMALL)
        {
            continue;
        }

        triad faceCoordSys(e0, e1, eN);
        faceCoordSys.normalize();

        // Construct the matrix to solve
        scalarSymmetricSquareMatrix T(3, 3, 0);
        scalarDiagonalMatrix Z(3, 0);

        // Least Squares
        for (label i = 0; i < 3; ++i)
        {
            scalar x = edgeVectors[i] & faceCoordSys[0];
            scalar y = edgeVectors[i] & faceCoordSys[1];

            T[0][0] += sqr(x);
            T[1][0] += x*y;
            T[1][1] += sqr(x) + sqr(y);
            T[2][1] += x*y;
            T[2][2] += sqr(y);

            scalar dndx = normalDifferences[i] & faceCoordSys[0];
            scalar dndy = normalDifferences[i] & faceCoordSys[1];

            Z[0] += dndx*x;
            Z[1] += dndx*y + dndy*x;
            Z[2] += dndy*y;
        }

        // Perform Cholesky decomposition and back substitution.
        // Decomposed matrix is in T and solution is in Z.
        LUsolve(T, Z);
        symmTensor2D secondFundamentalTensor(Z[0], Z[1], Z[2]);

        // Loop over the face points adding the contribution of the face
        // curvature to the points.
        forAll(f, fpI)
        {
            const label patchPointIndex = meshPointMap[f[fpI]];

            const triad& ptCoordSys = pointCoordSys[patchPointIndex];

            if (!ptCoordSys.set())
            {
                continue;
            }

            // Rotate faceCoordSys to ptCoordSys
            tensor rotTensor = rotationTensor(ptCoordSys[2], faceCoordSys[2]);
            triad rotatedFaceCoordSys = rotTensor & tensor(faceCoordSys);

            // Project the face curvature onto the point plane

            vector2D cmp1
            (
                ptCoordSys[0] & rotatedFaceCoordSys[0],
                ptCoordSys[0] & rotatedFaceCoordSys[1]
            );
            vector2D cmp2
            (
                ptCoordSys[1] & rotatedFaceCoordSys[0],
                ptCoordSys[1] & rotatedFaceCoordSys[1]
            );

            tensor2D projTensor
            (
                cmp1,
                cmp2
            );

            symmTensor2D projectedFundamentalTensor
            (
                projTensor.x() & (secondFundamentalTensor & projTensor.x()),
                projTensor.x() & (secondFundamentalTensor & projTensor.y()),
                projTensor.y() & (secondFundamentalTensor & projTensor.y())
            );

            // Calculate weight
            // @todo Voronoi area weighting
            scalar weight = calcVertexNormalWeight
            (
                f,
                meshPoints[patchPointIndex],
                f.normal(points),
                points
            );

            // Sum contribution of face to this point
            pointFundamentalTensors[patchPointIndex] +=
                weight*projectedFundamentalTensor;

            accumulatedWeights[patchPointIndex] += weight;
        }

        if (false)
        {
            Info<< "Points = " << points[f[0]] << " "
                << points[f[1]] << " "
                << points[f[2]] << endl;
            Info<< "edgeVecs = " << edgeVectors[0] << " "
                << edgeVectors[1] << " "
                << edgeVectors[2] << endl;
            Info<< "normDiff = " << normalDifferences[0] << " "
                << normalDifferences[1] << " "
                << normalDifferences[2] << endl;
            Info<< "faceCoordSys = " << faceCoordSys << endl;
            Info<< "T = " << T << endl;
            Info<< "Z = " << Z << endl;
        }
    }

    forAll(curvaturePointField, pI)
    {
        pointFundamentalTensors[pI] /= (accumulatedWeights[pI] + SMALL);

        vector2D principalCurvatures = eigenValues(pointFundamentalTensors[pI]);

        //scalar curvature =
        //    (principalCurvatures[0] + principalCurvatures[1])/2;
        scalar curvature = max
        (
            mag(principalCurvatures[0]),
            mag(principalCurvatures[1])
        );
        //scalar curvature = principalCurvatures[0]*principalCurvatures[1];

        curvaturePointField[meshPoints[pI]] = curvature;
    }

    return curvaturePointField;
}


bool edgesConnected(const edge& e1, const edge& e2)
{
    if
    (
        e1.start() == e2.start()
     || e1.start() == e2.end()
     || e1.end() == e2.start()
     || e1.end() == e2.end()
    )
    {
        return true;
    }

    return false;
}


scalar calcProximityOfFeaturePoints
(
    const List<pointIndexHit>& hitList,
    const scalar defaultCellSize
)
{
    scalar minDist = defaultCellSize;

    for
    (
        label hI1 = 0;
        hI1 < hitList.size() - 1;
        ++hI1
    )
    {
        const pointIndexHit& pHit1 = hitList[hI1];

        if (pHit1.hit())
        {
            for
            (
                label hI2 = hI1 + 1;
                hI2 < hitList.size();
                ++hI2
            )
            {
                const pointIndexHit& pHit2 = hitList[hI2];

                if (pHit2.hit())
                {
                    scalar curDist = mag(pHit1.hitPoint() - pHit2.hitPoint());

                    minDist = min(curDist, minDist);
                }
            }
        }
    }

    return minDist;
}


scalar calcProximityOfFeatureEdges
(
    const extendedFeatureEdgeMesh& efem,
    const List<pointIndexHit>& hitList,
    const scalar defaultCellSize
)
{
    scalar minDist = defaultCellSize;

    for
    (
        label hI1 = 0;
        hI1 < hitList.size() - 1;
        ++hI1
    )
    {
        const pointIndexHit& pHit1 = hitList[hI1];

        if (pHit1.hit())
        {
            const edge& e1 = efem.edges()[pHit1.index()];

            for
            (
                label hI2 = hI1 + 1;
                hI2 < hitList.size();
                ++hI2
            )
            {
                const pointIndexHit& pHit2 = hitList[hI2];

                if (pHit2.hit())
                {
                    const edge& e2 = efem.edges()[pHit2.index()];

                    // Don't refine if the edges are connected to each other
                    if (!edgesConnected(e1, e2))
                    {
                        scalar curDist =
                            mag(pHit1.hitPoint() - pHit2.hitPoint());

                        minDist = min(curDist, minDist);
                    }
                }
            }
        }
    }

    return minDist;
}


void dumpBox(const treeBoundBox& bb, const fileName& fName)
{
    OFstream str(fName);

    Info<< "Dumping bounding box " << bb << " as lines to obj file "
        << str.name() << endl;


    pointField boxPoints(bb.points());

    forAll(boxPoints, i)
    {
        meshTools::writeOBJ(str, boxPoints[i]);
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0]+1 <<  ' ' << e[1]+1 << nl;
    }
}


// Deletes all edges inside/outside bounding box from set.
void deleteBox
(
    const triSurface& surf,
    const treeBoundBox& bb,
    const bool removeInside,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgeI)
    {
        const point eMid = surf.edges()[edgeI].centre(surf.localPoints());

        if (removeInside ? bb.contains(eMid) : !bb.contains(eMid))
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


bool onLine(const point& p, const linePointRef& line)
{
    const point& a = line.start();
    const point& b = line.end();

    if
    (
        ( p.x() < min(a.x(), b.x()) || p.x() > max(a.x(), b.x()) )
     || ( p.y() < min(a.y(), b.y()) || p.y() > max(a.y(), b.y()) )
     || ( p.z() < min(a.z(), b.z()) || p.z() > max(a.z(), b.z()) )
    )
    {
        return false;
    }

    return true;
}


// Deletes all edges inside/outside bounding box from set.
void deleteEdges
(
    const triSurface& surf,
    const plane& cutPlane,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    const pointField& points = surf.points();
    const labelList& meshPoints = surf.meshPoints();

    forAll(edgeStat, edgeI)
    {
        const edge& e = surf.edges()[edgeI];
        const point& p0 = points[meshPoints[e.start()]];
        const point& p1 = points[meshPoints[e.end()]];
        const linePointRef line(p0, p1);

        // If edge does not intersect the plane, delete.
        scalar intersect = cutPlane.lineIntersect(line);

        point featPoint = intersect * (p1 - p0) + p0;

        if (!onLine(featPoint, line))
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


void drawHitProblem
(
    label fI,
    const triSurface& surf,
    const pointField& start,
    const pointField& faceCentres,
    const pointField& end,
    const List<pointIndexHit>& hitInfo
)
{
    Info<< nl << "# findLineAll did not hit its own face."
        << nl << "# fI " << fI
        << nl << "# start " << start[fI]
        << nl << "# f centre " << faceCentres[fI]
        << nl << "# end " << end[fI]
        << nl << "# hitInfo " << hitInfo
        << endl;

    meshTools::writeOBJ(Info, start[fI]);
    meshTools::writeOBJ(Info, faceCentres[fI]);
    meshTools::writeOBJ(Info, end[fI]);

    Info<< "l 1 2 3" << endl;

    meshTools::writeOBJ(Info, surf.points()[surf[fI][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][2]]);

    Info<< "f 4 5 6" << endl;

    forAll(hitInfo, hI)
    {
        label hFI = hitInfo[hI].index();

        meshTools::writeOBJ(Info, surf.points()[surf[hFI][0]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][1]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][2]]);

        Info<< "f "
            << 3*hI + 7 << " "
            << 3*hI + 8 << " "
            << 3*hI + 9
            << endl;
    }
}


// Unmark non-manifold edges if individual triangles are not features
void unmarkBaffles
(
    const triSurface& surf,
    const scalar includedAngle,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            label i0 = eFaces[0];
            //const labelledTri& f0 = surf[i0];
            const Foam::vector& n0 = surf.faceNormals()[i0];

            //Pout<< "edge:" << edgeI << " n0:" << n0 << endl;

            bool same = true;

            for (label i = 1; i < eFaces.size(); i++)
            {
                //const labelledTri& f = surf[i];
                const Foam::vector& n = surf.faceNormals()[eFaces[i]];

                //Pout<< "    mag(n&n0): " << mag(n&n0) << endl;

                if (mag(n&n0) < minCos)
                {
                    same = false;
                    break;
                }
            }

            if (same)
            {
                edgeStat[edgeI] = surfaceFeatures::NONE;
            }
        }
    }
}


//- Divide into multiple normal bins
//  - return REGION if != 2 normals
//  - return REGION if 2 normals that make feature angle
//  - otherwise return NONE and set normals,bins
surfaceFeatures::edgeStatus checkFlatRegionEdge
(
    const triSurface& surf,
    const scalar tol,
    const scalar includedAngle,
    const label edgeI
)
{
    const edge& e = surf.edges()[edgeI];
    const labelList& eFaces = surf.edgeFaces()[edgeI];

    // Bin according to normal

    DynamicList<Foam::vector> normals(2);
    DynamicList<labelList> bins(2);

    forAll(eFaces, eFaceI)
    {
        const Foam::vector& n = surf.faceNormals()[eFaces[eFaceI]];

        // Find the normal in normals
        label index = -1;
        forAll(normals, normalI)
        {
            if (mag(n&normals[normalI]) > (1-tol))
            {
                index = normalI;
                break;
            }
        }

        if (index != -1)
        {
            bins[index].append(eFaceI);
        }
        else if (normals.size() >= 2)
        {
            // Would be third normal. Mark as feature.
            //Pout<< "** at edge:" << surf.localPoints()[e[0]]
            //    << surf.localPoints()[e[1]]
            //    << " have normals:" << normals
            //    << " and " << n << endl;
            return surfaceFeatures::REGION;
        }
        else
        {
            normals.append(n);
            bins.append(labelList(1, eFaceI));
        }
    }


    // Check resulting number of bins
    if (bins.size() == 1)
    {
        // Note: should check here whether they are two sets of faces
        // that are planar or indeed 4 faces al coming together at an edge.
        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]]
        //    << " have single normal:" << normals[0]
        //    << endl;
        return surfaceFeatures::NONE;
    }
    else
    {
        // Two bins. Check if normals make an angle

        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]] << nl
        //    << "    normals:" << normals << nl
        //    << "    bins   :" << bins << nl
        //    << endl;

        if (includedAngle >= 0)
        {
            scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

            forAll(eFaces, i)
            {
                const Foam::vector& ni = surf.faceNormals()[eFaces[i]];
                for (label j=i+1; j<eFaces.size(); j++)
                {
                    const Foam::vector& nj = surf.faceNormals()[eFaces[j]];
                    if (mag(ni & nj) < minCos)
                    {
                        //Pout<< "have sharp feature between normal:" << ni
                        //    << " and " << nj << endl;

                        // Is feature. Keep as region or convert to
                        // feature angle? For now keep as region.
                        return surfaceFeatures::REGION;
                    }
                }
            }
        }


        // So now we have two normals bins but need to make sure both
        // bins have the same regions in it.

         // 1. store + or - region number depending
        //    on orientation of triangle in bins[0]
        const labelList& bin0 = bins[0];
        labelList regionAndNormal(bin0.size());
        forAll(bin0, i)
        {
            const labelledTri& t = surf.localFaces()[eFaces[bin0[i]]];
            int dir = t.edgeDirection(e);

            if (dir > 0)
            {
                regionAndNormal[i] = t.region()+1;
            }
            else if (dir == 0)
            {
                FatalErrorIn("problem.")
                    << exit(FatalError);
            }
            else
            {
                regionAndNormal[i] = -(t.region()+1);
            }
        }

        // 2. check against bin1
        const labelList& bin1 = bins[1];
        labelList regionAndNormal1(bin1.size());
        forAll(bin1, i)
        {
            const labelledTri& t = surf.localFaces()[eFaces[bin1[i]]];
            int dir = t.edgeDirection(e);

            label myRegionAndNormal;
            if (dir > 0)
            {
                myRegionAndNormal = t.region()+1;
            }
            else
            {
                myRegionAndNormal = -(t.region()+1);
            }

            regionAndNormal1[i] = myRegionAndNormal;

            label index = findIndex(regionAndNormal, -myRegionAndNormal);
            if (index == -1)
            {
                // Not found.
                //Pout<< "cannot find region " << myRegionAndNormal
                //    << " in regions " << regionAndNormal << endl;

                return surfaceFeatures::REGION;
            }
        }

        // Passed all checks, two normal bins with the same contents.
        //Pout<< "regionAndNormal:" << regionAndNormal << endl;
        //Pout<< "myRegionAndNormal:" << regionAndNormal1 << endl;

        return surfaceFeatures::NONE;
    }
}


void writeStats(const extendedFeatureEdgeMesh& fem, Ostream& os)
{
    os  << "    points : " << fem.points().size() << nl
        << "    of which" << nl
        << "        convex             : "
        << fem.concaveStart() << nl
        << "        concave            : "
        << (fem.mixedStart()-fem.concaveStart()) << nl
        << "        mixed              : "
        << (fem.nonFeatureStart()-fem.mixedStart()) << nl
        << "        non-feature        : "
        << (fem.points().size()-fem.nonFeatureStart()) << nl
        << "    edges  : " << fem.edges().size() << nl
        << "    of which" << nl
        << "        external edges     : "
        << fem.internalStart() << nl
        << "        internal edges     : "
        << (fem.flatStart()- fem.internalStart()) << nl
        << "        flat edges         : "
        << (fem.openStart()- fem.flatStart()) << nl
        << "        open edges         : "
        << (fem.multipleStart()- fem.openStart()) << nl
        << "        multiply connected : "
        << (fem.edges().size()- fem.multipleStart()) << endl;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();

#   include "addDictOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    const word dictName("surfaceFeatureExtractDict");
#   include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    const IOdictionary dict(dictIO);

    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            continue;
        }

        const dictionary& surfaceDict = iter().dict();

        if (!surfaceDict.found("extractionMethod"))
        {
            continue;
        }

        const word extractionMethod = surfaceDict.lookup("extractionMethod");

        const fileName surfFileName = iter().keyword();
        const fileName sFeatFileName = surfFileName.lessExt().name();

        Info<< "Surface            : " << surfFileName << nl << endl;

        const Switch writeVTK =
            surfaceDict.lookupOrDefault<Switch>("writeVTK", "off");
        const Switch writeObj =
            surfaceDict.lookupOrDefault<Switch>("writeObj", "off");

        const Switch curvature =
            surfaceDict.lookupOrDefault<Switch>("curvature", "off");
        const Switch featureProximity =
            surfaceDict.lookupOrDefault<Switch>("featureProximity", "off");
        const Switch closeness =
            surfaceDict.lookupOrDefault<Switch>("closeness", "off");


        Info<< nl << "Feature line extraction is only valid on closed manifold "
            << "surfaces." << endl;

        // Read
        // ~~~~

        triSurface surf(runTime.constantPath()/"triSurface"/surfFileName);

        Info<< "Statistics:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        faceList faces(surf.size());

        forAll(surf, fI)
        {
            faces[fI] = surf[fI].triFaceFace();
        }


        // Either construct features from surface & featureAngle or read set.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        autoPtr<surfaceFeatures> set;

        scalar includedAngle = 0.0;

        if (extractionMethod == "extractFromFile")
        {
            const dictionary& extractFromFileDict =
                surfaceDict.subDict("extractFromFileCoeffs");

            const fileName featureEdgeFile =
                extractFromFileDict.lookup("featureEdgeFile");

            const Switch geometricTestOnly =
                extractFromFileDict.lookupOrDefault<Switch>
                (
                    "geometricTestOnly",
                    "no"
                );

            edgeMesh eMesh(featureEdgeFile);

            // Sometimes duplicate edges are present. Remove them.
            eMesh.mergeEdges();

            Info<< nl << "Reading existing feature edges from file "
                << featureEdgeFile << nl
                << "Selecting edges purely based on geometric tests: "
                << geometricTestOnly.asText() << endl;

            set.set
            (
                new surfaceFeatures
                (
                    surf,
                    eMesh.points(),
                    eMesh.edges(),
                    1e-6,
                    geometricTestOnly
                )
            );
        }
        else if (extractionMethod == "extractFromSurface")
        {
            const dictionary& extractFromSurfaceDict =
                surfaceDict.subDict("extractFromSurfaceCoeffs");

            includedAngle =
                readScalar(extractFromSurfaceDict.lookup("includedAngle"));

            const Switch geometricTestOnly =
                extractFromSurfaceDict.lookupOrDefault<Switch>
                (
                    "geometricTestOnly",
                    "no"
                );

            Info<< nl
                << "Constructing feature set from included angle "
                << includedAngle << nl
                << "Selecting edges purely based on geometric tests: "
                << geometricTestOnly.asText() << endl;

            set.set
            (
                new surfaceFeatures
                (
                    surf,
                    includedAngle,
                    0,
                    0,
                    geometricTestOnly
                )
            );
        }
        else
        {
            FatalErrorIn(args.executable())
                << "No initial feature set. Provide either one"
                << " of extractFromFile (to read existing set)" << nl
                << " or extractFromSurface (to construct new set from angle)"
                << exit(FatalError);
        }


        // Trim set
        // ~~~~~~~~

        if (surfaceDict.isDict("trimFeatures"))
        {
            dictionary trimDict = surfaceDict.subDict("trimFeatures");

            scalar minLen =
                trimDict.lookupOrAddDefault<scalar>("minLen", -GREAT);

            label minElem = trimDict.lookupOrAddDefault<label>("minElem", 0);

            // Trim away small groups of features
            if (minElem > 0 || minLen > 0)
            {
                Info<< "Removing features of length < "
                    << minLen << endl;
                Info<< "Removing features with number of edges < "
                    << minElem << endl;

                set().trimFeatures(minLen, minElem, includedAngle);
            }
        }


        // Subset
        // ~~~~~~

        // Convert to marked edges, points
        List<surfaceFeatures::edgeStatus> edgeStat(set().toStatus());

        if (surfaceDict.isDict("subsetFeatures"))
        {
            const dictionary& subsetDict = surfaceDict.subDict
            (
                "subsetFeatures"
            );

            if (subsetDict.found("insideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("insideBox")());

                Info<< "Removing all edges outside bb " << bb << endl;
                dumpBox(bb, "subsetBox.obj");

                deleteBox(surf, bb, false, edgeStat);
            }
            else if (subsetDict.found("outsideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("outsideBox")());

                Info<< "Removing all edges inside bb " << bb << endl;
                dumpBox(bb, "deleteBox.obj");

                deleteBox(surf, bb, true, edgeStat);
            }

            const Switch nonManifoldEdges =
                subsetDict.lookupOrDefault<Switch>("nonManifoldEdges", "yes");

            if (!nonManifoldEdges)
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                forAll(edgeStat, edgeI)
                {
                    const labelList& eFaces = surf.edgeFaces()[edgeI];

                    if
                    (
                        eFaces.size() > 2
                     && edgeStat[edgeI] == surfaceFeatures::REGION
                     && (eFaces.size() % 2) == 0
                    )
                    {
                        edgeStat[edgeI] = checkFlatRegionEdge
                        (
                            surf,
                            1e-5,   //tol,
                            includedAngle,
                            edgeI
                        );
                    }
                }
            }

            const Switch openEdges =
                subsetDict.lookupOrDefault<Switch>("openEdges", "yes");

            if (!openEdges)
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                forAll(edgeStat, edgeI)
                {
                    if (surf.edgeFaces()[edgeI].size() == 1)
                    {
                        edgeStat[edgeI] = surfaceFeatures::NONE;
                    }
                }
            }

            if (subsetDict.found("plane"))
            {
                plane cutPlane(subsetDict.lookup("plane")());

                deleteEdges(surf, cutPlane, edgeStat);

                Info<< "Only edges that intersect the plane with normal "
                    << cutPlane.normal()
                    << " and base point " << cutPlane.refPoint()
                    << " will be included as feature edges."<< endl;
            }
        }


        surfaceFeatures newSet(surf);
        newSet.setFromStatus(edgeStat, includedAngle);

        Info<< nl
            << "Initial feature set:" << nl
            << "    feature points : " << newSet.featurePoints().size() << nl
            << "    feature edges  : " << newSet.featureEdges().size() << nl
            << "    of which" << nl
            << "        region edges   : " << newSet.nRegionEdges() << nl
            << "        external edges : " << newSet.nExternalEdges() << nl
            << "        internal edges : " << newSet.nInternalEdges() << nl
            << endl;

        //if (writeObj)
        //{
        //    newSet.writeObj("final");
        //}

        boolList surfBaffleRegions(surf.patches().size(), false);

        wordList surfBaffleNames;
        surfaceDict.readIfPresent("baffles", surfBaffleNames);

        forAll(surf.patches(), pI)
        {
            const word& name = surf.patches()[pI].name();

            if (findIndex(surfBaffleNames, name) != -1)
            {
                Info<< "Adding baffle region " << name << endl;
                surfBaffleRegions[pI] = true;
            }
        }

        // Extracting and writing a extendedFeatureEdgeMesh
        extendedFeatureEdgeMesh feMesh
        (
            newSet,
            runTime,
            sFeatFileName + ".extendedFeatureEdgeMesh",
            surfBaffleRegions
        );


        if (surfaceDict.isDict("addFeatures"))
        {
            const word addFeName = surfaceDict.subDict("addFeatures")["name"];
            Info<< "Adding (without merging) features from " << addFeName
                << nl << endl;

            extendedFeatureEdgeMesh addFeMesh
            (
                IOobject
                (
                    addFeName,
                    runTime.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            Info<< "Read " << addFeMesh.name() << nl;
            writeStats(addFeMesh, Info);

            feMesh.add(addFeMesh);
        }


        Info<< nl
            << "Final feature set:" << nl;
        writeStats(feMesh, Info);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.objectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj(feMesh.path()/surfFileName.lessExt().name());
        }

        feMesh.write();

        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfFileName.lessExt().name() + ".eMesh",   // name
                runTime.constant(),                         // instance
                "triSurface",
                runTime,                                    // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectPath() << endl;

        bfeMesh.regIOobject::write();

    // Find close features

    // // Dummy trim operation to mark features
    // labelList featureEdgeIndexing = newSet.trimFeatures(-GREAT, 0);

    // scalarField surfacePtFeatureIndex(surf.points().size(), -1);

    // forAll(newSet.featureEdges(), eI)
    // {
    //     const edge& e = surf.edges()[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.start()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.end()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];
    // }

    // if (writeVTK)
    // {
    //     vtkSurfaceWriter().write
    //     (
    //         runTime.constant()/"triSurface",    // outputDir
    //         sFeatFileName,                      // surfaceName
    //         surf.points(),
    //         faces,
    //         "surfacePtFeatureIndex",            // fieldName
    //         surfacePtFeatureIndex,
    //         true,                               // isNodeValues
    //         true                                // verbose
    //     );
    // }

    // Random rndGen(343267);

    // treeBoundBox surfBB
    // (
    //     treeBoundBox(searchSurf.bounds()).extend(rndGen, 1e-4)
    // );

    // surfBB.min() -= Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    // surfBB.max() += Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // indexedOctree<treeDataEdge> ftEdTree
    // (
    //     treeDataEdge
    //     (
    //         false,
    //         surf.edges(),
    //         surf.localPoints(),
    //         newSet.featureEdges()
    //     ),
    //     surfBB,
    //     8,      // maxLevel
    //     10,     // leafsize
    //     3.0     // duplicity
    // );

    // labelList nearPoints = ftEdTree.findBox
    // (
    //     treeBoundBox
    //     (
    //         sPt - featureSearchSpan*Foam::vector::one,
    //         sPt + featureSearchSpan*Foam::vector::one
    //     )
    // );

        if (closeness)
        {
            Info<< nl << "Extracting internal and external closeness of "
                << "surface." << endl;


            triSurfaceMesh searchSurf
            (
                IOobject
                (
                    sFeatFileName + ".closeness",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf
            );


            // Internal and external closeness

            // Prepare start and end points for intersection tests

            const vectorField& normals = searchSurf.faceNormals();

            scalar span = searchSurf.bounds().mag();

            scalar externalAngleTolerance = 10;
            scalar externalToleranceCosAngle =
                Foam::cos
                (
                    degToRad(180 - externalAngleTolerance)
                );

            scalar internalAngleTolerance = 45;
            scalar internalToleranceCosAngle =
                Foam::cos
                (
                    degToRad(180 - internalAngleTolerance)
                );

            Info<< "externalToleranceCosAngle: " << externalToleranceCosAngle
                << nl
                << "internalToleranceCosAngle: " << internalToleranceCosAngle
                << endl;

            // Info<< "span " << span << endl;

            pointField start(searchSurf.faceCentres() - span*normals);
            pointField end(searchSurf.faceCentres() + span*normals);
            const pointField& faceCentres = searchSurf.faceCentres();

            List<List<pointIndexHit> > allHitInfo;

            // Find all intersections (in order)
            searchSurf.findLineAll(start, end, allHitInfo);

            scalarField internalCloseness(start.size(), GREAT);
            scalarField externalCloseness(start.size(), GREAT);

            forAll(allHitInfo, fI)
            {
                const List<pointIndexHit>& hitInfo = allHitInfo[fI];

                if (hitInfo.size() < 1)
                {
                    drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

                    // FatalErrorIn(args.executable())
                    //     << "findLineAll did not hit its own face."
                    //     << exit(FatalError);
                }
                else if (hitInfo.size() == 1)
                {
                    if (!hitInfo[0].hit())
                    {
                        // FatalErrorIn(args.executable())
                        //     << "findLineAll did not hit any face."
                        //     << exit(FatalError);
                    }
                    else if (hitInfo[0].index() != fI)
                    {
                        drawHitProblem
                        (
                            fI,
                            surf,
                            start,
                            faceCentres,
                            end,
                            hitInfo
                        );

                        // FatalErrorIn(args.executable())
                        //     << "findLineAll did not hit its own face."
                        //     << exit(FatalError);
                    }
                }
                else
                {
                    label ownHitI = -1;

                    forAll(hitInfo, hI)
                    {
                        // Find the hit on the triangle that launched the ray

                        if (hitInfo[hI].index() == fI)
                        {
                            ownHitI = hI;

                            break;
                        }
                    }

                    if (ownHitI < 0)
                    {
                        drawHitProblem
                        (
                            fI,
                            surf,
                            start,
                            faceCentres,
                            end,
                            hitInfo
                        );

                        // FatalErrorIn(args.executable())
                        //     << "findLineAll did not hit its own face."
                        //     << exit(FatalError);
                    }
                    else if (ownHitI == 0)
                    {
                        // There are no internal hits, the first hit is the
                        // closest external hit

                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI + 1].index()]
                            )
                          < externalToleranceCosAngle
                        )
                        {
                            externalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI + 1].hitPoint()
                                );
                        }
                    }
                    else if (ownHitI == hitInfo.size() - 1)
                    {
                        // There are no external hits, the last but one hit is
                        // the closest internal hit

                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI - 1].index()]
                            )
                          < internalToleranceCosAngle
                        )
                        {
                            internalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI - 1].hitPoint()
                                );
                        }
                    }
                    else
                    {
                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI + 1].index()]
                            )
                          < externalToleranceCosAngle
                        )
                        {
                            externalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI + 1].hitPoint()
                                );
                        }

                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI - 1].index()]
                            )
                          < internalToleranceCosAngle
                        )
                        {
                            internalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI - 1].hitPoint()
                                );
                        }
                    }
                }
            }

            triSurfaceScalarField internalClosenessField
            (
                IOobject
                (
                    sFeatFileName + ".internalCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                internalCloseness
            );

            internalClosenessField.write();

            triSurfaceScalarField externalClosenessField
            (
                IOobject
                (
                    sFeatFileName + ".externalCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                externalCloseness
            );

            externalClosenessField.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "internalCloseness",                // fieldName
                    internalCloseness,
                    false,                              // isNodeValues
                    true                                // verbose
                );

                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "externalCloseness",                // fieldName
                    externalCloseness,
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }


        if (curvature)
        {
            Info<< nl << "Extracting curvature of surface at the points."
                << endl;

            vectorField pointNormals = calcVertexNormals(surf);
            triadField pointCoordSys = calcVertexCoordSys(surf, pointNormals);

            triSurfacePointScalarField k = calcCurvature
            (
                sFeatFileName,
                runTime,
                surf,
                pointNormals,
                pointCoordSys
            );

            k.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "curvature",                        // fieldName
                    k,
                    true,                               // isNodeValues
                    true                                // verbose
                );
            }
        }


        if (featureProximity)
        {
            Info<< nl << "Extracting proximity of close feature points and "
                << "edges to the surface" << endl;

            const scalar searchDistance =
                readScalar(surfaceDict.lookup("maxFeatureProximity"));

            scalarField featureProximity(surf.size(), searchDistance);

            forAll(surf, fI)
            {
                const triPointRef& tri = surf[fI].tri(surf.points());
                const point& triCentre = tri.circumCentre();

                const scalar radiusSqr = min
                (
                    sqr(4*tri.circumRadius()),
                    sqr(searchDistance)
                );

                List<pointIndexHit> hitList;

                feMesh.allNearestFeatureEdges(triCentre, radiusSqr, hitList);

                featureProximity[fI] =
                    calcProximityOfFeatureEdges
                    (
                        feMesh,
                        hitList,
                        featureProximity[fI]
                    );

                feMesh.allNearestFeaturePoints(triCentre, radiusSqr, hitList);

                featureProximity[fI] =
                    calcProximityOfFeaturePoints
                    (
                        hitList,
                        featureProximity[fI]
                    );
            }

            triSurfaceScalarField featureProximityField
            (
                IOobject
                (
                    sFeatFileName + ".featureProximity",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                featureProximity
            );

            featureProximityField.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "featureProximity",                 // fieldName
                    featureProximity,
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }

        Info<< endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
