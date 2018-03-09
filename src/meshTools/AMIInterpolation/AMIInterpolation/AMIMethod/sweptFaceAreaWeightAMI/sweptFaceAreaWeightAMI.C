/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "sweptFaceAreaWeightAMI.H"
#include "cut.H"
#include "linearEqn.H"
#include "quadraticEqn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
const Foam::scalar Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::
minCutRatio_ = 10*small;

template<class SourcePatch, class TargetPatch>
const Foam::scalar Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::
maxDot_ = - cos(degToRad(89.9));


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
template<unsigned Size>
void Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::writeCutTrisVTK
(
    const cutTriList<Size>& tris,
    const word& name
) const
{
    OFstream obj(typeName + "_" + name + ".vtk");

    obj << "# vtk DataFile Version 2.0" << endl
        << "triSurface" << endl
        << "ASCII" << endl
        << "DATASET POLYDATA" << endl
        << "POINTS " << 3*tris.size() << " float" << endl;

    for (label i = 0; i < tris.size(); ++ i)
    {
        for (label j = 0; j < 3; ++ j)
        {
            const vector& x = tris[i][j];
            obj << x.x() << ' ' << x.y() << ' ' << x.z() << endl;
        }
    }

    obj << "POLYGONS " << tris.size() << ' ' << 4*tris.size() << endl;

    for (label i = 0; i < tris.size(); ++ i)
    {
        obj << 3 << ' ' << 3*i << ' ' << 3*i + 1 << ' ' << 3*i + 2 << endl;
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::writeFaceOBJ
(
    const face& f,
    const pointField& ps,
    const string& name
) const
{
    OFstream obj(typeName + "_" + name + ".obj");

    for (label i = 0; i < f.size(); ++ i)
    {
        const vector& x = ps[f[i]];
        obj << "v " << x.x() << ' ' << x.y() << ' ' << x.z() << endl;
    }

    obj << 'f';
    for (label i = 0; i < f.size(); ++ i)
    {
        obj << ' ' << i + 1;
    }
    obj << endl;
}


template<class SourcePatch, class TargetPatch>
void Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::writeProjectionOBJ
(
    const label srcN,
    const FixedList<point, 4>& srcTri,
    const FixedList<point, 4>& srcNrm
) const
{
    scalar l = 0;
    for (label i = 0; i < srcN - 1; ++ i)
    {
        l = max(l, mag(srcTri[i] - srcTri[i + 1]));
    }

    const label nu = 20, nv = 20;
    const scalar u0 = 0, u1 = 1, v0 = -l, v1 = l;

    OFstream obj(typeName + "_projection.obj");

    for (label i = 0; i < srcN; ++ i)
    {
        const point& p0 = srcTri[i], p1 = srcTri[(i + 1) % srcN];
        const vector& n0 = srcNrm[i], n1 = srcNrm[(i + 1) % srcN];

        for (label iu = 0; iu <= nu; ++ iu)
        {
            const scalar u = u0 + (u1 - u0)*scalar(iu)/nu;
            for (label iv = 0; iv <= nv; ++ iv)
            {
                const scalar v = v0 + (v1 - v0)*scalar(iv)/nv;
                const vector x = p0 + (p1 - p0)*u + (n0 + (n1 - n0)*u)*v;
                obj << "v " << x.x() << ' ' << x.y() << ' ' << x.z() << endl;
            }
        }
    }

    for (label i = 0; i < srcN; ++ i)
    {
        for (label iu = 0; iu < nu; ++ iu)
        {
            for (label iv = 0; iv < nv; ++ iv)
            {
                obj << "f "
                    << i*(nu + 1)*(nv + 1) + (nv + 1)*iu + iv + 1 << ' '
                    << i*(nu + 1)*(nv + 1) + (nv + 1)*iu + iv + 2 << ' '
                    << i*(nu + 1)*(nv + 1) + (nv + 1)*(iu + 1) + iv + 2 << ' '
                    << i*(nu + 1)*(nv + 1) + (nv + 1)*(iu + 1) + iv + 1<< endl;
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label
Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::getSourceProjection
(
    FixedList<point, 4>& srcTri,
    FixedList<point, 4>& srcNrm,
    const FixedList<point, 3>& tgtTri
) const
{
    // The target normal
    const vector np = plane(tgtTri[0], tgtTri[1], tgtTri[2]).normal();

    // Dot products between the projection normals and the target plane
    FixedList<scalar, 4> dots;
    for (label i = 0; i < 3; ++ i)
    {
        dots[i] = srcNrm[i] & np;
    }
    dots[3] = 0;

    // Search though for the number of reversed directions and an index for
    // both a reverse and a forward direction
    label nR = 0, iR = -1, iF = -1;
    for (label i = 0; i < 3; ++ i)
    {
        if (dots[i] > maxDot_)
        {
            ++ nR;
            iR = i;
        }
        else
        {
            iF = i;
        }
    }

    // If all the normals hit in the forward or reverse direction then return
    // the existing projection or no projection respectively
    if (nR == 0 || nR == 3)
    {
        return 3 - nR;
    }

    // One normal hits in the reverse direction
    if (nR == 1)
    {
        /*
                       o j1
                      / \
                     /   \
                    /     \
                   /       o jR
                  /       / .
                 /       /   .
                /       /     .
               o - > - o . . . . (old iR)
              i0      iR
        */

        // Make a gap by duplicating the reverse the point
        for (label j = 2; j >= iR; -- j)
        {
            srcTri[j + 1] = srcTri[j];
            srcNrm[j + 1] = srcNrm[j];
            dots[j + 1] = dots[j];
        }

        const label jR = iR + 1;
        const label i0 = (iR + 3) % 4, j1 = (jR + 1) % 4;

        const scalar wi = (dots[iR] - maxDot_)/(dots[iR] - dots[i0]);
        const scalar wj = (dots[jR] - maxDot_)/(dots[jR] - dots[j1]);

        srcTri[iR] = srcTri[iR] + (srcTri[i0] - srcTri[iR])*wi;
        srcTri[jR] = srcTri[jR] + (srcTri[j1] - srcTri[jR])*wj;

        srcNrm[iR] = srcNrm[iR] + (srcNrm[i0] - srcNrm[iR])*wi;
        srcNrm[jR] = srcNrm[jR] + (srcNrm[j1] - srcNrm[jR])*wj;

        return 4;
    }

    // Two normals hit in the reverse direction
    if (nR == 2)
    {
        /*
                       . (old jR)
                      . .
                  jR o   .
                    / \   .
                   /   \   .
                  /     \   .
                 /       \   .
                /         \   .
               o - - > - - o . . (old iR)
              iF          iR
        */

        const label iR = (iF + 1) % 3, jR = (iF + 2) % 3;
        const label i0 = (iR + 2) % 3, j1 = (jR + 1) % 3;

        const scalar wi = (dots[iR] - maxDot_)/(dots[iR] - dots[i0]);
        const scalar wj = (dots[jR] - maxDot_)/(dots[jR] - dots[j1]);

        srcTri[iR] = srcTri[iR] + (srcTri[i0] - srcTri[iR])*wi;
        srcTri[jR] = srcTri[jR] + (srcTri[j1] - srcTri[jR])*wj;

        srcNrm[iR] = srcNrm[iR] + (srcNrm[i0] - srcNrm[iR])*wi;
        srcNrm[jR] = srcNrm[jR] + (srcNrm[j1] - srcNrm[jR])*wj;

        return 3;
    }

    // This cannot happen
    return -1;
}


template<class SourcePatch, class TargetPatch>
Foam::plane Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::getCutPlane
(
    const point& p0,
    const point& p1,
    const vector& n0,
    const vector& n1,
    const FixedList<point, 3>& tgtTri
) const
{
    // The target plane
    const plane tgtPln(tgtTri[0], tgtTri[1], tgtTri[2]);
    const point& pp = tgtPln.refPoint();
    const vector& np = tgtPln.normal();

    // Calculate the bounding intersection points. These are the locations at
    // which the bounding lines of the projected surface intersect with the
    // target plane.
    const scalar v0 = ((pp - p0) & np)/(n0 & np);
    const scalar v1 = ((pp - p1) & np)/(n1 & np);
    Pair<point> cutPnts(p0 + v0*n0, p1 + v1*n1);
    Pair<vector> cutNrms((p1 - p0 + n1*v0) ^ n0, (p1 - p0 - n0*v1) ^ n1);

    // Calculate edge intersections with the surface
    for (label i = 0; i < 3; ++ i)
    {
        // The current target edge
        const point& l0 = tgtTri[i], l1 = tgtTri[(i + 1)%3];

        // Form the quadratic in the surface parameter u
        const vector k = l0 - p0;
        const vector a = l0 - l1, b = p1 - p0, c = n0, d = n1 - n0;
        const vector ka = k ^ a, ba = b ^ a, ca = c ^ a, da = d ^ a;
        const scalar A = d & ba;
        const scalar B = (c & ba) - (d & ka);
        const scalar C = - c & ka;

        // Solve and extract the other parameters
        const Roots<2> us = quadraticEqn(A, B, C).roots();
        Pair<scalar> vs, ws;
        Pair<vector> ns;
        Pair<bool> cuts(false, false);
        for (label j = 0; j < 2; ++ j)
        {
            if (us.type(j) == roots::real)
            {
                const vector den = ca + da*us[j];
                if (magSqr(den) > vSmall)
                {
                    const vector vNum = ka - ba*us[j];
                    const vector wNum = (- k + b*us[j]) ^ (c + d*us[j]);
                    vs[j] = (vNum & den)/magSqr(den);
                    ws[j] = (wNum & den)/magSqr(den);
                    ns[j] = (b ^ (c + d*us[j])) + (n1 ^ n0)*vs[j];
                    const bool cutu = 0 < us[j] && us[j] < 1;
                    const bool cutw = 0 < ws[j] && ws[j] < 1;
                    cuts[j] = cutu && cutw;
                }
            }
        }

        // If we have just one intersection in bounds then store it in the
        // result list based on its direction
        if (cuts[0] != cuts[1])
        {
            const label j = cuts[0] ? 0 : 1;
            const scalar na = ns[j] & a;
            cutPnts[na < 0] = l0 + (l1 - l0)*ws[j];
            cutNrms[na < 0] = ns[j];
        }
    }

    // Generate and return the cut plane. If the cut points are not coincident
    // then form a plane normal by crossing the displacement between the points
    // by the target plane normal. If the points are coincident then use the
    // projected surface normal evaluated at the first cut point.
    const vector cutDelta = cutPnts[1] - cutPnts[0];
    const bool coincident = mag(cutDelta) < minCutRatio_*mag(p1 - p0);
    return plane(cutPnts[0], coincident ? cutNrms[0] : np ^ cutDelta);
};


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::scalar Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::interArea
(
    const label srcFacei,
    const label tgtFacei
) const
{
    const label debugTgtFace =
        (- 1 - debug) % this->tgtPatch_.size() == tgtFacei
      ? (- 1 - debug) / this->tgtPatch_.size() + 1 : 0;
    const bool debugPrint = debug > 0 || debugTgtFace > 0;
    const bool debugWrite = debug > 1 || debugTgtFace > 1;

    if (debugPrint)
    {
        Info<< "Inter area between source face #" << srcFacei
            << " and target face #" << tgtFacei << (debugWrite ? "\n" : ": ");
    }

    // Patch data
    const pointField& srcPoints = this->srcPatch_.localPoints();
    const pointField& tgtPoints = this->tgtPatch_.localPoints();
    const vectorField& srcPointNormals = this->srcPatch_.pointNormals();

    // Faces
    const face& srcFace = this->srcPatch_.localFaces()[srcFacei];
    const face& tgtFace = this->tgtPatch_.localFaces()[tgtFacei];

    // Write out the faces
    if (debugWrite)
    {
        writeFaceOBJ(srcFace, srcPoints, "source");
        writeFaceOBJ(tgtFace, tgtPoints, "target");
    }
    if (debugTgtFace)
    {
        writeFaceOBJ(tgtFace, tgtPoints, "target" + name(tgtFacei));
    }

    // Triangulate the faces
    const faceAreaIntersect::triangulationMode triMode = this->triMode_;
    faceList srcFaceTris, tgtFaceTris;
    faceAreaIntersect::triangulate(srcFace, srcPoints, triMode, srcFaceTris);
    faceAreaIntersect::triangulate(tgtFace, tgtPoints, triMode, tgtFaceTris);

    // Area sum
    scalar areaMag = Zero;

    // Loop the target triangles
    forAllConstIter(faceList, tgtFaceTris, tgtIter)
    {
        const FixedList<point, 3>
            tgtTri =
            {
                tgtPoints[(*tgtIter)[0]],
                tgtPoints[(*tgtIter)[this->reverseTarget_ ? 2 : 1]],
                tgtPoints[(*tgtIter)[this->reverseTarget_ ? 1 : 2]]
            };

        // Loop the source triangles
        forAllConstIter(faceList, srcFaceTris, srcIter)
        {
            FixedList<point, 4>
                srcTri =
                {
                    srcPoints[(*srcIter)[0]],
                    srcPoints[(*srcIter)[1]],
                    srcPoints[(*srcIter)[2]],
                    vector::zero
                };
            FixedList<point, 4>
                srcNrm =
                {
                    srcPointNormals[(*srcIter)[0]],
                    srcPointNormals[(*srcIter)[1]],
                    srcPointNormals[(*srcIter)[2]],
                    vector::zero
                };

            // Get the source projection modified for any reverse intersections
            const label srcN = getSourceProjection(srcTri, srcNrm, tgtTri);
            if (srcN <= 0)
            {
                continue;
            }

            // Write the source projection
            if (debugWrite)
            {
                writeProjectionOBJ(srcN, srcTri, srcNrm);
            }

            // Create the initial cut triangle list
            cutTriList<8> cutTris;
            cutTris.append(tgtTri);

            // Write the initial triangle
            if (debugWrite)
            {
                writeCutTrisVTK(cutTris, "tris0");
            }

            // Do all but one of the cuts
            for
            (
                label i = 0;
                i < srcN - 1 && (debugWrite || cutTris.size());
                ++ i
            )
            {
                // Do the cut
                const plane cutPlane =
                    getCutPlane
                    (
                        srcTri[i],
                        srcTri[i+1],
                        srcNrm[i],
                        srcNrm[i+1],
                        tgtTri
                    );

                cutTriList<8> cutTrisTmp;

                for (label j = 0; j < cutTris.size(); ++j)
                {
                    triCut
                    (
                        cutTris[j],
                        cutPlane,
                        cut::noOp(),
                        cut::appendOp<cutTriList<8>>(cutTrisTmp)
                    );
                }
                Swap(cutTris, cutTrisTmp);

                // Write the triangles resulting from the cut
                if (debugWrite)
                {
                    writeCutTrisVTK(cutTris, "tris" + name(i + 1));
                }
            }

            // Do the last cut
            const plane cutPlane =
                getCutPlane
                (
                    srcTri[srcN - 1],
                    srcTri[0],
                    srcNrm[srcN - 1],
                    srcNrm[0],
                    tgtTri
                );
            cutTriList<8> cutTrisTmp;
            for (label i = 0; i < cutTris.size(); ++ i)
            {
                // Sum the area of the cut triangle
                areaMag +=
                    mag
                    (
                        triCut
                        (
                            cutTris[i],
                            cutPlane,
                            cut::noOp(),
                            cut::areaOp()
                        )
                    );
                // Store the cut triangle if it is needed for output
                if (debugWrite || debugTgtFace)
                {
                    triCut
                    (
                        cutTris[i],
                        cutPlane,
                        cut::noOp(),
                        cut::appendOp<cutTriList<8>>(cutTrisTmp)
                    );
                }
            }
            Swap(cutTris, cutTrisTmp);

            // Write the triangles resulting from the cuts
            if (debugTgtFace && cutTris.size())
            {
                static label writeCutTrisVTKIndex = 0;
                writeCutTrisVTK
                (
                    cutTris,
                    "target" + name(tgtFacei) + "_"
                  + "tris" + name(writeCutTrisVTKIndex ++)
                );
            }
            if (debugWrite)
            {
                writeCutTrisVTK(cutTris, "tris" + name(srcN));

                Info << "view triangles then press ENTER to continue ...";
                getchar();
            }
        }
    }

    // Print the difference between this inter-area and that obtained by the
    // standard algorithm built around faceAreaIntersect
    if (debugPrint)
    {
        const scalar standardAreaMag =
            faceAreaWeightAMI<SourcePatch, TargetPatch>::interArea
            (
                srcFacei,
                tgtFacei
            );

        Info<<"standard=" << standardAreaMag << ", swept=" << areaMag
            << endl;
    }

    return areaMag;
}


template<class SourcePatch, class TargetPatch>
Foam::scalar
Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::minWeight() const
{
    return small;
}


template<class SourcePatch, class TargetPatch>
Foam::scalar
Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::maxWalkAngle() const
{
    return degToRad(180);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::sweptFaceAreaWeightAMI<SourcePatch, TargetPatch>::~sweptFaceAreaWeightAMI
(
)
{}


// ************************************************************************* //
