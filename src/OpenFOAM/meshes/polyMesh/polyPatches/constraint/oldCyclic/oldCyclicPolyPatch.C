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

\*---------------------------------------------------------------------------*/

#include "oldCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "patchZones.H"
#include "matchPoints.H"
#include "Time.H"
#include "transformList.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oldCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, oldCyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, oldCyclicPolyPatch, dictionary);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointField Foam::oldCyclicPolyPatch::calcFaceCentres
(
    const UList<face>& faces,
    const pointField& points
)
{
    pointField ctrs(faces.size());

    forAll(faces, facei)
    {
        ctrs[facei] = faces[facei].centre(points);
    }

    return ctrs;
}


Foam::pointField Foam::oldCyclicPolyPatch::getAnchorPoints
(
    const UList<face>& faces,
    const pointField& points
)
{
    pointField anchors(faces.size());

    forAll(faces, facei)
    {
        anchors[facei] = points[faces[facei][0]];
    }

    return anchors;
}


Foam::label Foam::oldCyclicPolyPatch::findMaxArea
(
    const pointField& points,
    const faceList& faces
)
{
    label maxI = -1;
    scalar maxAreaSqr = -great;

    forAll(faces, facei)
    {
        scalar areaSqr = magSqr(faces[facei].area(points));

        if (areaSqr > maxAreaSqr)
        {
            maxAreaSqr = areaSqr;
            maxI = facei;
        }
    }
    return maxI;
}


bool Foam::oldCyclicPolyPatch::getGeometricHalves
(
    const primitivePatch& pp,
    labelList& half0ToPatch,
    labelList& half1ToPatch
) const
{
    // Get geometric zones of patch by looking at normals.
    // Method 1: any edge with sharpish angle is edge between two halves.
    //           (this will handle e.g. wedge geometries).
    //           Also two fully disconnected regions will be handled this way.
    // Method 2: sort faces into two halves based on face normal.

    // Calculate normals
    const vectorField& faceNormals = pp.faceNormals();

    // Find edges with sharp angles.
    boolList regionEdge(pp.nEdges(), false);

    const labelListList& edgeFaces = pp.edgeFaces();

    label nRegionEdges = 0;

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        // Check manifold edges for sharp angle.
        // (Non-manifold already handled by patchZones)
        if (eFaces.size() == 2)
        {
            if ((faceNormals[eFaces[0]] & faceNormals[eFaces[1]])< featureCos_)
            {
                regionEdge[edgeI] = true;

                nRegionEdges++;
            }
        }
    }


    // For every face determine zone it is connected to (without crossing
    // any regionEdge)
    patchZones ppZones(pp, regionEdge);

    if (debug)
    {
        Pout<< "oldCyclicPolyPatch::getGeometricHalves : "
            << "Found " << nRegionEdges << " edges on patch " << name()
            << " where the cos of the angle between two connected faces"
            << " was less than " << featureCos_ << nl
            << "Patch divided by these and by single sides edges into "
            << ppZones.nZones() << " parts." << endl;


        // Dumping zones to obj files.

        labelList nZoneFaces(ppZones.nZones());

        for (label zoneI = 0; zoneI < ppZones.nZones(); zoneI++)
        {
            OFstream stream
            (
                boundaryMesh().mesh().time().path()
               /name()+"_zone_"+Foam::name(zoneI)+".obj"
            );
            Pout<< "oldCyclicPolyPatch::getGeometricHalves : Writing zone "
                << zoneI << " face centres to OBJ file " << stream.name()
                << endl;

            labelList zoneFaces(findIndices(ppZones, zoneI));

            forAll(zoneFaces, i)
            {
                writeOBJ(stream, pp[zoneFaces[i]].centre(pp.points()));
            }

            nZoneFaces[zoneI] = zoneFaces.size();
        }
    }


    if (ppZones.nZones() == 2)
    {
        half0ToPatch = findIndices(ppZones, 0);
        half1ToPatch = findIndices(ppZones, 1);
    }
    else
    {
        if (debug)
        {
            Pout<< "oldCyclicPolyPatch::getGeometricHalves :"
                << " falling back to face-normal comparison" << endl;
        }
        label n0Faces = 0;
        half0ToPatch.setSize(pp.size());

        label n1Faces = 0;
        half1ToPatch.setSize(pp.size());

        // Compare to face 0 normal.
        forAll(faceNormals, facei)
        {
            if ((faceNormals[facei] & faceNormals[0]) > 0)
            {
                half0ToPatch[n0Faces++] = facei;
            }
            else
            {
                half1ToPatch[n1Faces++] = facei;
            }
        }
        half0ToPatch.setSize(n0Faces);
        half1ToPatch.setSize(n1Faces);

        if (debug)
        {
            Pout<< "oldCyclicPolyPatch::getGeometricHalves :"
                << " Number of faces per zone:("
                << n0Faces << ' ' << n1Faces << ')' << endl;
        }
    }

    if (half0ToPatch.size() != half1ToPatch.size())
    {
        fileName casePath(boundaryMesh().mesh().time().path());

        // Dump halves
        {
            fileName nm0(casePath/name()+"_half0_faces.obj");
            Pout<< "oldCyclicPolyPatch::getGeometricHalves : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, UIndirectList<face>(pp, half0ToPatch)(), pp.points());

            fileName nm1(casePath/name()+"_half1_faces.obj");
            Pout<< "oldCyclicPolyPatch::getGeometricHalves : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, UIndirectList<face>(pp, half1ToPatch)(), pp.points());
        }

        // Dump face centres
        {
            OFstream str0(casePath/name()+"_half0.obj");
            Pout<< "oldCyclicPolyPatch::getGeometricHalves : Writing half0"
                << " face centres to OBJ file " << str0.name() << endl;

            forAll(half0ToPatch, i)
            {
                writeOBJ(str0, pp[half0ToPatch[i]].centre(pp.points()));
            }

            OFstream str1(casePath/name()+"_half1.obj");
            Pout<< "oldCyclicPolyPatch::getGeometricHalves : Writing half1"
                << " face centres to OBJ file " << str1.name() << endl;
            forAll(half1ToPatch, i)
            {
                writeOBJ(str1, pp[half1ToPatch[i]].centre(pp.points()));
            }
        }

        SeriousErrorInFunction
            << "Patch " << name() << " gets decomposed in two zones of"
            << "inequal size: " << half0ToPatch.size()
            << " and " << half1ToPatch.size() << endl
            << "This means that the patch is either not two separate regions"
            << " or one region where the angle between the different regions"
            << " is not sufficiently sharp." << endl
            << "Please adapt the featureCos setting." << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::oldCyclicPolyPatch::getCentresAndAnchors
(
    const primitivePatch& pp,
    const faceList& half0Faces,
    const faceList& half1Faces,

    pointField& ppPoints,
    pointField& half0Ctrs,
    pointField& half1Ctrs,
    pointField& anchors0,
    scalarField& tols
) const
{
    // Get geometric data on both halves.
    half0Ctrs = calcFaceCentres(half0Faces, pp.points());
    anchors0 = getAnchorPoints(half0Faces, pp.points());
    half1Ctrs = calcFaceCentres(half1Faces, pp.points());

    switch (transform())
    {
        case ROTATIONAL:
        {
            label face0 = getConsistentRotationFace(half0Ctrs);
            label face1 = getConsistentRotationFace(half1Ctrs);

            vector n0 = ((half0Ctrs[face0] - rotationCentre_) ^ rotationAxis_);
            vector n1 = ((half1Ctrs[face1] - rotationCentre_) ^ -rotationAxis_);
            n0 /= mag(n0) + vSmall;
            n1 /= mag(n1) + vSmall;

            if (debug)
            {
                Pout<< "oldCyclicPolyPatch::getCentresAndAnchors :"
                    << " Specified rotation :"
                    << " n0:" << n0 << " n1:" << n1 << endl;
            }

            // Rotation (around origin)
            const tensor reverseT(rotationTensor(n0, -n1));

            // Rotation
            forAll(half0Ctrs, facei)
            {
                half0Ctrs[facei] = Foam::transform(reverseT, half0Ctrs[facei]);
                anchors0[facei] = Foam::transform(reverseT, anchors0[facei]);
            }

            ppPoints = Foam::transform(reverseT, pp.points());

            break;
        }
        //- Problem: usually specified translation is not accurate enough
        //- To get proper match so keep automatic determination over here.
        //case TRANSLATIONAL:
        //{
        //    // Transform 0 points.
        //
        //    if (debug)
        //    {
        //        Pout<< "oldCyclicPolyPatch::getCentresAndAnchors :"
        //            << "Specified translation : " << separationVector_
        //            << endl;
        //    }
        //
        //    half0Ctrs += separationVector_;
        //    anchors0 += separationVector_;
        //    break;
        //}
        default:
        {
            // Assumes that cyclic is planar. This is also the initial
            // condition for patches without faces.

            // Determine the face with max area on both halves. These
            // two faces are used to determine the transformation tensors
            const label max0I = findMaxArea(pp.points(), half0Faces);
            const vector n0 = half0Faces[max0I].normal(pp.points());

            const label max1I = findMaxArea(pp.points(), half1Faces);
            const vector n1 = half1Faces[max1I].normal(pp.points());

            if (mag(n0 & n1) < 1-matchTolerance())
            {
                if (debug)
                {
                    Pout<< "oldCyclicPolyPatch::getCentresAndAnchors :"
                        << " Detected rotation :"
                        << " n0:" << n0 << " n1:" << n1 << endl;
                }

                // Rotation (around origin)
                const tensor reverseT(rotationTensor(n0, -n1));

                // Rotation
                forAll(half0Ctrs, facei)
                {
                    half0Ctrs[facei] = Foam::transform
                    (
                        reverseT,
                        half0Ctrs[facei]
                    );
                    anchors0[facei] = Foam::transform
                    (
                        reverseT,
                        anchors0[facei]
                    );
                }
                ppPoints = Foam::transform(reverseT, pp.points());
            }
            else
            {
                // Parallel translation. Get average of all used points.

                primitiveFacePatch half0(half0Faces, pp.points());
                const pointField& half0Pts = half0.localPoints();
                const point ctr0(sum(half0Pts)/half0Pts.size());

                primitiveFacePatch half1(half1Faces, pp.points());
                const pointField& half1Pts = half1.localPoints();
                const point ctr1(sum(half1Pts)/half1Pts.size());

                if (debug)
                {
                    Pout<< "oldCyclicPolyPatch::getCentresAndAnchors :"
                        << " Detected translation :"
                        << " n0:" << n0 << " n1:" << n1
                        << " ctr0:" << ctr0 << " ctr1:" << ctr1 << endl;
                }

                half0Ctrs += ctr1 - ctr0;
                anchors0 += ctr1 - ctr0;
                ppPoints = pp.points() + ctr1 - ctr0;
            }
            break;
        }
    }


    // Calculate typical distance per face
    tols = matchTolerance()*calcFaceTol(half1Faces, pp.points(), half1Ctrs);
}


bool Foam::oldCyclicPolyPatch::matchAnchors
(
    const bool report,
    const primitivePatch& pp,
    const labelList& half0ToPatch,
    const pointField& anchors0,

    const labelList& half1ToPatch,
    const faceList& half1Faces,
    const labelList& from1To0,

    const scalarField& tols,

    labelList& faceMap,
    labelList& rotation
) const
{
    // Set faceMap such that half0 faces get first and corresponding half1
    // faces last.

    forAll(half0ToPatch, half0Facei)
    {
        // Label in original patch
        label patchFacei = half0ToPatch[half0Facei];

        faceMap[patchFacei] = half0Facei;

        // No rotation
        rotation[patchFacei] = 0;
    }

    bool fullMatch = true;

    forAll(from1To0, half1Facei)
    {
        label patchFacei = half1ToPatch[half1Facei];

        // This face has to match the corresponding one on half0.
        label half0Facei = from1To0[half1Facei];

        label newFacei = half0Facei + pp.size()/2;

        faceMap[patchFacei] = newFacei;

        // Rotate patchFacei such that its f[0] aligns with that of
        // the corresponding face
        // (which after shuffling will be at position half0Facei)

        const point& wantedAnchor = anchors0[half0Facei];

        rotation[newFacei] = getRotation
        (
            pp.points(),
            half1Faces[half1Facei],
            wantedAnchor,
            tols[half1Facei]
        );

        if (rotation[newFacei] == -1)
        {
            fullMatch = false;

            if (report)
            {
                const face& f = half1Faces[half1Facei];
                SeriousErrorInFunction
                    << "Patch:" << name() << " : "
                    << "Cannot find point on face " << f
                    << " with vertices:"
                    << UIndirectList<point>(pp.points(), f)()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of cyclic patch " << name()
                    << endl
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;
            }
        }
    }
    return fullMatch;
}


Foam::label Foam::oldCyclicPolyPatch::getConsistentRotationFace
(
    const pointField& faceCentres
) const
{
    const scalarField magRadSqr
    (
        magSqr((faceCentres - rotationCentre_) ^ rotationAxis_)
    );
    scalarField axisLen
    (
        (faceCentres - rotationCentre_) & rotationAxis_
    );
    axisLen = axisLen - min(axisLen);
    const scalarField magLenSqr
    (
        magRadSqr + axisLen*axisLen
    );

    label rotFace = -1;
    scalar maxMagLenSqr = -great;
    scalar maxMagRadSqr = -great;
    forAll(faceCentres, i)
    {
        if (magLenSqr[i] >= maxMagLenSqr)
        {
            if (magRadSqr[i] > maxMagRadSqr)
            {
                rotFace = i;
                maxMagLenSqr = magLenSqr[i];
                maxMagRadSqr = magRadSqr[i];
            }
        }
    }

    if (debug)
    {
        Info<< "getConsistentRotationFace(const pointField&)" << nl
            << "    rotFace = " << rotFace << nl
            << "    point =  " << faceCentres[rotFace] << endl;
    }

    return rotFace;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::oldCyclicPolyPatch::oldCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    featureCos_(0.9),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    separationVector_(Zero)
{}


Foam::oldCyclicPolyPatch::oldCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    featureCos_(0.9),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    separationVector_(Zero)
{
    if (dict.found("neighbourPatch"))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Found \"neighbourPatch\" entry when reading cyclic patch "
            << name << endl
            << "Is this mesh already with split cyclics?" << endl
            << "If so run a newer version that supports it"
            << ", if not comment out the \"neighbourPatch\" entry and re-run"
            << exit(FatalIOError);
    }

    dict.readIfPresent("featureCos", featureCos_);

    switch (transform())
    {
        case ROTATIONAL:
        {
            dict.lookup("rotationAxis") >> rotationAxis_;
            dict.lookup("rotationCentre") >> rotationCentre_;
            break;
        }
        case TRANSLATIONAL:
        {
            dict.lookup("separationVector") >> separationVector_;
            break;
        }
        default:
        {
            // no additional info required
        }
    }
}


Foam::oldCyclicPolyPatch::oldCyclicPolyPatch
(
    const oldCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    featureCos_(pp.featureCos_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_)
{}


Foam::oldCyclicPolyPatch::oldCyclicPolyPatch
(
    const oldCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    featureCos_(pp.featureCos_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oldCyclicPolyPatch::~oldCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oldCyclicPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::oldCyclicPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{}


void Foam::oldCyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{}


void Foam::oldCyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::oldCyclicPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
}


void Foam::oldCyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::oldCyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


void Foam::oldCyclicPolyPatch::initOrder
(
    PstreamBuffers&,
    const primitivePatch& pp
) const
{}


//  Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool Foam::oldCyclicPolyPatch::order
(
    PstreamBuffers&,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (pp.empty())
    {
        // No faces, nothing to change.
        return false;
    }

    if (pp.size()&1)
    {
        FatalErrorInFunction
            << "Size of cyclic " << name() << " should be a multiple of 2"
            << ". It is " << pp.size() << abort(FatalError);
    }

    label halfSize = pp.size()/2;

    // Supplied primitivePatch already with new points.
    // Cyclics are limited to one transformation tensor
    // currently anyway (i.e. straight plane) so should not be too big a
    // problem.


    // Indices of faces on half0
    labelList half0ToPatch;
    // Indices of faces on half1
    labelList half1ToPatch;


    // 1. Test if already correctly oriented by starting from trivial ordering.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    half0ToPatch = identity(halfSize);
    half1ToPatch = half0ToPatch + halfSize;

    // Get faces
    faceList half0Faces(UIndirectList<face>(pp, half0ToPatch));
    faceList half1Faces(UIndirectList<face>(pp, half1ToPatch));

    // Get geometric quantities
    pointField half0Ctrs, half1Ctrs, anchors0, ppPoints;
    scalarField tols;
    getCentresAndAnchors
    (
        pp,
        half0Faces,
        half1Faces,

        ppPoints,
        half0Ctrs,
        half1Ctrs,
        anchors0,
        tols
    );

    // Geometric match of face centre vectors
    labelList from1To0;
    bool matchedAll = matchPoints
    (
        half1Ctrs,
        half0Ctrs,
        tols,
        false,
        from1To0
    );

    if (debug)
    {
        Pout<< "oldCyclicPolyPatch::order : test if already ordered:"
            << matchedAll << endl;

        // Dump halves
        fileName nm0("match1_"+name()+"_half0_faces.obj");
        Pout<< "oldCyclicPolyPatch::order : Writing half0"
            << " faces to OBJ file " << nm0 << endl;
        writeOBJ(nm0, half0Faces, ppPoints);

        fileName nm1("match1_"+name()+"_half1_faces.obj");
        Pout<< "oldCyclicPolyPatch::order : Writing half1"
            << " faces to OBJ file " << nm1 << endl;
        writeOBJ(nm1, half1Faces, pp.points());

        OFstream ccStr
        (
            boundaryMesh().mesh().time().path()
           /"match1_"+ name() + "_faceCentres.obj"
        );
        Pout<< "oldCyclicPolyPatch::order : "
            << "Dumping currently found cyclic match as lines between"
            << " corresponding face centres to file " << ccStr.name()
            << endl;

        // Recalculate untransformed face centres
        //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
        label vertI = 0;

        forAll(half1Ctrs, i)
        {
            //if (from1To0[i] != -1)
            {
                // Write edge between c1 and c0
                //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                //const point& c0 = half0Ctrs[from1To0[i]];
                const point& c0 = half0Ctrs[i];
                const point& c1 = half1Ctrs[i];
                writeOBJ(ccStr, c0, c1, vertI);
            }
        }
    }


    // 2. Ordered in pairs (so 0,1 coupled and 2,3 etc.)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!matchedAll)
    {
        label facei = 0;
        for (label i = 0; i < halfSize; i++)
        {
            half0ToPatch[i] = facei++;
            half1ToPatch[i] = facei++;
        }

        // And redo all matching
        half0Faces = UIndirectList<face>(pp, half0ToPatch);
        half1Faces = UIndirectList<face>(pp, half1ToPatch);

        getCentresAndAnchors
        (
            pp,
            half0Faces,
            half1Faces,

            ppPoints,
            half0Ctrs,
            half1Ctrs,
            anchors0,
            tols
        );

        // Geometric match of face centre vectors
        matchedAll = matchPoints
        (
            half1Ctrs,
            half0Ctrs,
            tols,
            false,
            from1To0
        );

        if (debug)
        {
            Pout<< "oldCyclicPolyPatch::order : test if pairwise ordered:"
                << matchedAll << endl;
            // Dump halves
            fileName nm0("match2_"+name()+"_half0_faces.obj");
            Pout<< "oldCyclicPolyPatch::order : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0Faces, ppPoints);

            fileName nm1("match2_"+name()+"_half1_faces.obj");
            Pout<< "oldCyclicPolyPatch::order : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1Faces, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /"match2_"+name()+"_faceCentres.obj"
            );
            Pout<< "oldCyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            label vertI = 0;

            forAll(half1Ctrs, i)
            {
                if (from1To0[i] != -1)
                {
                    // Write edge between c1 and c0
                    const point& c0 = half0Ctrs[from1To0[i]];
                    const point& c1 = half1Ctrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }
    }


    // 3. Baffles(coincident faces) converted into cyclics (e.g. jump)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!matchedAll)
    {
        label baffleI = 0;

        forAll(pp, facei)
        {
            const face& f = pp.localFaces()[facei];
            const labelList& pFaces = pp.pointFaces()[f[0]];

            label matchedFacei = -1;

            forAll(pFaces, i)
            {
                label otherFacei = pFaces[i];

                if (otherFacei > facei)
                {
                    const face& otherF = pp.localFaces()[otherFacei];

                    // Note: might pick up two similar oriented faces
                    //       (but that is illegal anyway)
                    if (f == otherF)
                    {
                        matchedFacei = otherFacei;
                        break;
                    }
                }
            }

            if (matchedFacei != -1)
            {
                half0ToPatch[baffleI] = facei;
                half1ToPatch[baffleI] = matchedFacei;
                baffleI++;
            }
        }

        if (baffleI == halfSize)
        {
            // And redo all matching
            half0Faces = UIndirectList<face>(pp, half0ToPatch);
            half1Faces = UIndirectList<face>(pp, half1ToPatch);

            getCentresAndAnchors
            (
                pp,
                half0Faces,
                half1Faces,

                ppPoints,
                half0Ctrs,
                half1Ctrs,
                anchors0,
                tols
            );

            // Geometric match of face centre vectors
            matchedAll = matchPoints
            (
                half1Ctrs,
                half0Ctrs,
                tols,
                false,
                from1To0
            );

            if (debug)
            {
                Pout<< "oldCyclicPolyPatch::order : test if baffles:"
                    << matchedAll << endl;
                // Dump halves
                fileName nm0("match3_"+name()+"_half0_faces.obj");
                Pout<< "oldCyclicPolyPatch::order : Writing half0"
                    << " faces to OBJ file " << nm0 << endl;
                writeOBJ(nm0, half0Faces, ppPoints);

                fileName nm1("match3_"+name()+"_half1_faces.obj");
                Pout<< "oldCyclicPolyPatch::order : Writing half1"
                    << " faces to OBJ file " << nm1 << endl;
                writeOBJ(nm1, half1Faces, pp.points());

                OFstream ccStr
                (
                    boundaryMesh().mesh().time().path()
                   /"match3_"+ name() + "_faceCentres.obj"
                );
                Pout<< "oldCyclicPolyPatch::order : "
                    << "Dumping currently found cyclic match as lines between"
                    << " corresponding face centres to file " << ccStr.name()
                    << endl;

                // Recalculate untransformed face centres
                label vertI = 0;

                forAll(half1Ctrs, i)
                {
                    if (from1To0[i] != -1)
                    {
                        // Write edge between c1 and c0
                        const point& c0 = half0Ctrs[from1To0[i]];
                        const point& c1 = half1Ctrs[i];
                        writeOBJ(ccStr, c0, c1, vertI);
                    }
                }
            }
        }
    }


    // 4. Automatic geometric ordering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!matchedAll)
    {
        // Split faces according to feature angle or topology
        bool okSplit = getGeometricHalves(pp, half0ToPatch, half1ToPatch);

        if (!okSplit)
        {
            // Did not split into two equal parts.
            return false;
        }

        // And redo all matching
        half0Faces = UIndirectList<face>(pp, half0ToPatch);
        half1Faces = UIndirectList<face>(pp, half1ToPatch);

        getCentresAndAnchors
        (
            pp,
            half0Faces,
            half1Faces,

            ppPoints,
            half0Ctrs,
            half1Ctrs,
            anchors0,
            tols
        );

        // Geometric match of face centre vectors
        matchedAll = matchPoints
        (
            half1Ctrs,
            half0Ctrs,
            tols,
            false,
            from1To0
        );

        if (debug)
        {
            Pout<< "oldCyclicPolyPatch::order : automatic ordering result:"
                << matchedAll << endl;
            // Dump halves
            fileName nm0("match4_"+name()+"_half0_faces.obj");
            Pout<< "oldCyclicPolyPatch::order : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0Faces, ppPoints);

            fileName nm1("match4_"+name()+"_half1_faces.obj");
            Pout<< "oldCyclicPolyPatch::order : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1Faces, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /"match4_"+ name() + "_faceCentres.obj"
            );
            Pout<< "oldCyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            label vertI = 0;

            forAll(half1Ctrs, i)
            {
                if (from1To0[i] != -1)
                {
                    // Write edge between c1 and c0
                    const point& c0 = half0Ctrs[from1To0[i]];
                    const point& c1 = half1Ctrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }
    }


    if (!matchedAll || debug)
    {
        // Dump halves
        fileName nm0(name()+"_half0_faces.obj");
        Pout<< "oldCyclicPolyPatch::order : Writing half0"
            << " faces to OBJ file " << nm0 << endl;
        writeOBJ(nm0, half0Faces, pp.points());

        fileName nm1(name()+"_half1_faces.obj");
        Pout<< "oldCyclicPolyPatch::order : Writing half1"
            << " faces to OBJ file " << nm1 << endl;
        writeOBJ(nm1, half1Faces, pp.points());

        OFstream ccStr
        (
            boundaryMesh().mesh().time().path()
           /name() + "_faceCentres.obj"
        );
        Pout<< "oldCyclicPolyPatch::order : "
            << "Dumping currently found cyclic match as lines between"
            << " corresponding face centres to file " << ccStr.name()
            << endl;

        // Recalculate untransformed face centres
        //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
        label vertI = 0;

        forAll(half1Ctrs, i)
        {
            if (from1To0[i] != -1)
            {
                // Write edge between c1 and c0
                //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                const point& c0 = half0Ctrs[from1To0[i]];
                const point& c1 = half1Ctrs[i];
                writeOBJ(ccStr, c0, c1, vertI);
            }
        }
    }


    if (!matchedAll)
    {
        SeriousErrorInFunction
            << "Patch:" << name() << " : "
            << "Cannot match vectors to faces on both sides of patch" << endl
            << "    Perhaps your faces do not match?"
            << " The obj files written contain the current match." << endl
            << "    Continuing with incorrect face ordering from now on!"
            << endl;

            return false;
    }


    // Set faceMap such that half0 faces get first and corresponding half1
    // faces last.
    matchAnchors
    (
        true,                   // report if anchor matching error
        pp,
        half0ToPatch,
        anchors0,
        half1ToPatch,
        half1Faces,
        from1To0,
        tols,
        faceMap,
        rotation
    );

    // Return false if no change necessary, true otherwise.

    forAll(faceMap, facei)
    {
        if (faceMap[facei] != facei || rotation[facei] != 0)
        {
            return true;
        }
    }

    return false;
}


void Foam::oldCyclicPolyPatch::write(Ostream& os) const
{
    // Replacement of polyPatch::write to write 'cyclic' instead of type():
    os.writeKeyword("type") << cyclicPolyPatch::typeName
        << token::END_STATEMENT << nl;
    patchIdentifier::write(os);
    os.writeKeyword("nFaces") << size() << token::END_STATEMENT << nl;
    os.writeKeyword("startFace") << start() << token::END_STATEMENT << nl;


    os.writeKeyword("featureCos") << featureCos_ << token::END_STATEMENT << nl;
    switch (transform())
    {
        case ROTATIONAL:
        {
            os.writeKeyword("rotationAxis") << rotationAxis_
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationCentre") << rotationCentre_
                << token::END_STATEMENT << nl;
            break;
        }
        case TRANSLATIONAL:
        {
            os.writeKeyword("separationVector") << separationVector_
                << token::END_STATEMENT << nl;
            break;
        }
        default:
        {
            // no additional info to write
        }
    }

    WarningInFunction
        << "Please run foamUpgradeCyclics to convert these old-style"
        << " cyclics into two separate cyclics patches."
        << endl;
}


// ************************************************************************* //
