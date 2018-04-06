/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    surfaceFeatureExtract

Description
    Extracts and writes surface features to file. All but the basic feature
    extraction is WIP.

    Curvature calculation is an implementation of the algorithm from:

      "Estimating Curvatures and their Derivatives on Triangle Meshes"
      by S. Rusinkiewicz

\*---------------------------------------------------------------------------*/

#include "surfaceFeatureExtract.H"
#include "Time.H"
#include "meshTools.H"
#include "tensor2D.H"
#include "symmTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::scalar Foam::internalAngleTolerance(80);
const Foam::scalar Foam::internalToleranceCosAngle
(
    cos(degToRad(180 - internalAngleTolerance))
);

const Foam::scalar Foam::externalAngleTolerance(10);
const Foam::scalar Foam::externalToleranceCosAngle
(
    cos(degToRad(180 - externalAngleTolerance))
);


Foam::point Foam::randomPointInPlane(const plane& p)
{
    // Perturb base point
    const point& refPt = p.refPoint();

    // ax + by + cz + d = 0
    const FixedList<scalar, 4>& planeCoeffs = p.planeCoeffs();

    const scalar perturbX = refPt.x() + 1e-3;
    const scalar perturbY = refPt.y() + 1e-3;
    const scalar perturbZ = refPt.z() + 1e-3;

    if (mag(planeCoeffs[2]) < small)
    {
        if (mag(planeCoeffs[1]) < small)
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


Foam::triadField Foam::calcVertexCoordSys
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

        if (mag(normal) < small)
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


Foam::triSurfacePointScalarField Foam::calcCurvature
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

    forAll(surf, fi)
    {
        const triFace& f = surf[fi];
        const edgeList fEdges = f.edges();

        // Calculate the edge vectors and the normal differences
        vectorField edgeVectors(f.size(), Zero);
        vectorField normalDifferences(f.size(), Zero);

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
        const vector eN = f.area(points);
        const vector e1 = (e0 ^ eN);

        if (magSqr(eN) < rootVSmall)
        {
            continue;
        }

        triad faceCoordSys(e0, e1, eN);
        faceCoordSys.normalize();

        // Construct the matrix to solve
        scalarSymmetricSquareMatrix T(3, 0);
        scalarDiagonalMatrix Z(3, 0);

        // Least Squares
        for (label i = 0; i < 3; ++i)
        {
            scalar x = edgeVectors[i] & faceCoordSys[0];
            scalar y = edgeVectors[i] & faceCoordSys[1];

            T(0, 0) += sqr(x);
            T(1, 0) += x*y;
            T(1, 1) += sqr(x) + sqr(y);
            T(2, 1) += x*y;
            T(2, 2) += sqr(y);

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
            // TODO: Voronoi area weighting
            const scalar weight = surf.pointNormalWeight
            (
                f,
                meshPoints[patchPointIndex],
                f.area(points),
                points
            );

            // Sum contribution of face to this point
            pointFundamentalTensors[patchPointIndex] +=
                weight*projectedFundamentalTensor;

            accumulatedWeights[patchPointIndex] += weight;
        }
    }

    forAll(curvaturePointField, pI)
    {
        pointFundamentalTensors[pI] /= (accumulatedWeights[pI] + small);

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


bool Foam::edgesConnected(const edge& e1, const edge& e2)
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


Foam::scalar Foam::calcProximityOfFeaturePoints
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


Foam::scalar Foam::calcProximityOfFeatureEdges
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


void Foam::dumpBox(const treeBoundBox& bb, const fileName& fName)
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


void Foam::deleteBox
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


bool Foam::onLine(const point& p, const linePointRef& line)
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


void Foam::deleteEdges
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


void Foam::drawHitProblem
(
    const label fi,
    const triSurface& surf,
    const point& start,
    const point& p,
    const point& end,
    const List<pointIndexHit>& hitInfo
)
{
    Info<< nl << "# findLineAll did not hit its own face."
        << nl << "# fi " << fi
        << nl << "# start " << start
        << nl << "# point " << p
        << nl << "# end " << end
        << nl << "# hitInfo " << hitInfo
        << endl;

    meshTools::writeOBJ(Info, start);
    meshTools::writeOBJ(Info, p);
    meshTools::writeOBJ(Info, end);

    Info<< "l 1 2 3" << endl;

    meshTools::writeOBJ(Info, surf.points()[surf[fi][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fi][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fi][2]]);

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


void Foam::unmarkBaffles
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


Foam::surfaceFeatures::edgeStatus Foam::checkFlatRegionEdge
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

    forAll(eFaces, eFacei)
    {
        const Foam::vector& n = surf.faceNormals()[eFaces[eFacei]];

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
            bins[index].append(eFacei);
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
            bins.append(labelList(1, eFacei));
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
                FatalErrorInFunction
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


void Foam::extractCloseness
(
    const fileName &sFeatFileName,
    const Time& runTime,
    const triSurface &surf,
    const bool writeVTK
);


void Foam::extractPointCloseness
(
    const fileName &sFeatFileName,
    const Time& runTime,
    const triSurface &surf,
    const bool writeVTK
);


void Foam::writeStats(const extendedFeatureEdgeMesh& fem, Ostream& os)
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


// ************************************************************************* //
