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

\*---------------------------------------------------------------------------*/

#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "matchPoints.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "Time.H"
#include "transformList.H"
#include "PstreamBuffers.H"
#include "ConstCirculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const transformType transform,
    const word& patchType
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const transformType transform,
    const word& patchType
)
:
    coupledPolyPatch
    (
        newName(myProcNo, neighbProcNo),
        size,
        start,
        index,
        bm,
        patchType,
        transform
    ),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::~processorPolyPatch()
{
    neighbPointsPtr_.clear();
    neighbEdgesPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::processorPolyPatch::newName
(
    const label myProcNo,
    const label neighbProcNo
)
{
    return
        "procBoundary"
      + Foam::name(myProcNo)
      + "to"
      + Foam::name(neighbProcNo);
}


void Foam::processorPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << faceCentres()
            << faceAreas()
            << faceCellCentres();
    }
}


void Foam::processorPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        {
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> neighbFaceCentres_
                >> neighbFaceAreas_
                >> neighbFaceCellCentres_;
        }

        // My normals
        vectorField faceNormals(size());

        // Neighbour normals
        vectorField nbrFaceNormals(neighbFaceAreas_.size());

        // Face match tolerances
        scalarField tols = calcFaceTol(*this, points(), faceCentres());

        // Calculate normals from areas and check
        forAll(faceNormals, facei)
        {
            scalar magSf = mag(faceAreas()[facei]);
            scalar nbrMagSf = mag(neighbFaceAreas_[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            // For small face area calculation the results of the area
            // calculation have been found to only be accurate to ~1e-20
            if (magSf < small || nbrMagSf < small)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check.
                faceNormals[facei] = point(1, 0, 0);
                nbrFaceNormals[facei] = -faceNormals[facei];
                tols[facei] = great;
            }
            else if (mag(magSf - nbrMagSf) > matchTolerance()*sqr(tols[facei]))
            {
                fileName nm
                (
                    boundaryMesh().mesh().time().path()
                   /name()+"_faces.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry : Writing my "
                    << size()
                    << " faces to OBJ file " << nm << endl;

                writeOBJ(nm, *this, points());

                OFstream ccStr
                (
                    boundaryMesh().mesh().time().path()
                    /name() + "_faceCentresConnections.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry :"
                    << " Dumping cell centre lines between"
                    << " corresponding face centres to OBJ file" << ccStr.name()
                    << endl;

                label vertI = 0;

                forAll(faceCentres(), facej)
                {
                    const point& c0 = neighbFaceCentres_[facej];
                    const point& c1 = faceCentres()[facej];

                    writeOBJ(ccStr, c0, c1, vertI);
                }

                FatalErrorInFunction
                    << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf - nbrMagSf)/avSf
                    << "% -- possible face ordering problem." << endl
                    << "patch:" << name()
                    << " my area:" << magSf
                    << " neighbour area:" << nbrMagSf
                    << " matching tolerance:"
                    << matchTolerance()*sqr(tols[facei])
                    << endl
                    << "Mesh face:" << start()+facei
                    << " vertices:"
                    << UIndirectList<point>(points(), operator[](facei))()
                    << endl
                    << "If you are certain your matching is correct"
                    << " you can increase the 'matchTolerance' setting"
                    << " in the patch dictionary in the boundary file."
                    << endl
                    << "Rerun with processor debug flag set for"
                    << " more information." << exit(FatalError);
            }
            else
            {
                faceNormals[facei] = faceAreas()[facei]/magSf;
                nbrFaceNormals[facei] = neighbFaceAreas_[facei]/nbrMagSf;
            }
        }

        calcTransformTensors
        (
            faceCentres(),
            neighbFaceCentres_,
            faceNormals,
            nbrFaceNormals,
            matchTolerance()*tols,
            matchTolerance(),
            transform()
        );
    }
}


void Foam::processorPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    processorPolyPatch::initGeometry(pBufs);
}


void Foam::processorPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField&
)
{
    processorPolyPatch::calcGeometry(pBufs);
}


void Foam::processorPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);

    if (Pstream::parRun())
    {
        // Express all points as patch face and index in face.
        labelList pointFace(nPoints());
        labelList pointIndex(nPoints());

        for (label patchPointi = 0; patchPointi < nPoints(); patchPointi++)
        {
            label facei = pointFaces()[patchPointi][0];

            pointFace[patchPointi] = facei;

            const face& f = localFaces()[facei];

            pointIndex[patchPointi] = findIndex(f, patchPointi);
        }

        // Express all edges as patch face and index in face.
        labelList edgeFace(nEdges());
        labelList edgeIndex(nEdges());

        for (label patchEdgeI = 0; patchEdgeI < nEdges(); patchEdgeI++)
        {
            label facei = edgeFaces()[patchEdgeI][0];

            edgeFace[patchEdgeI] = facei;

            const labelList& fEdges = faceEdges()[facei];

            edgeIndex[patchEdgeI] = findIndex(fEdges, patchEdgeI);
        }

        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << pointFace
            << pointIndex
            << edgeFace
            << edgeIndex;
    }
}


void Foam::processorPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    polyPatch::updateMesh(pBufs);

    neighbPointsPtr_.clear();
    neighbEdgesPtr_.clear();

    if (Pstream::parRun())
    {
        labelList nbrPointFace;
        labelList nbrPointIndex;
        labelList nbrEdgeFace;
        labelList nbrEdgeIndex;

        {
            // !Note: there is one situation where the opposite points and
            // do not exactly match and that is while we are distributing
            // meshes and multiple parts come together from different
            // processors. This can temporarily create the situation that
            // we have points which will be merged out later but we still
            // need the face connectivity to be correct.
            // So: cannot check here on same points and edges.

            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> nbrPointFace
                >> nbrPointIndex
                >> nbrEdgeFace
                >> nbrEdgeIndex;
        }

        // Convert neighbour faces and indices into face back into
        // my edges and points.

        // Convert points.
        // ~~~~~~~~~~~~~~~

        neighbPointsPtr_.reset(new labelList(nPoints(), -1));
        labelList& neighbPoints = neighbPointsPtr_();

        forAll(nbrPointFace, nbrPointi)
        {
            // Find face and index in face on this side.
            const face& f = localFaces()[nbrPointFace[nbrPointi]];

            label index = (f.size() - nbrPointIndex[nbrPointi]) % f.size();
            label patchPointi = f[index];

            if (neighbPoints[patchPointi] == -1)
            {
                // First reference of point
                neighbPoints[patchPointi] = nbrPointi;
            }
            else if (neighbPoints[patchPointi] >= 0)
            {
                // Point already visited. Mark as duplicate.
                neighbPoints[patchPointi] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(neighbPoints, patchPointi)
        {
            if (neighbPoints[patchPointi] == -2)
            {
                neighbPoints[patchPointi] = -1;
            }
        }

        // Convert edges.
        // ~~~~~~~~~~~~~~

        neighbEdgesPtr_.reset(new labelList(nEdges(), -1));
        labelList& neighbEdges = neighbEdgesPtr_();

        forAll(nbrEdgeFace, nbrEdgeI)
        {
            // Find face and index in face on this side.
            const labelList& f = faceEdges()[nbrEdgeFace[nbrEdgeI]];
            label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
            label patchEdgeI = f[index];

            if (neighbEdges[patchEdgeI] == -1)
            {
                // First reference of edge
                neighbEdges[patchEdgeI] = nbrEdgeI;
            }
            else if (neighbEdges[patchEdgeI] >= 0)
            {
                // Edge already visited. Mark as duplicate.
                neighbEdges[patchEdgeI] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(neighbEdges, patchEdgeI)
        {
            if (neighbEdges[patchEdgeI] == -2)
            {
                neighbEdges[patchEdgeI] = -1;
            }
        }

        // Remove any addressing used for shared points/edges calculation
        // since mostly not needed.
        primitivePatch::clearOut();
    }
}


const Foam::labelList& Foam::processorPolyPatch::neighbPoints() const
{
    if (!neighbPointsPtr_.valid())
    {
        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return neighbPointsPtr_();
}


const Foam::labelList& Foam::processorPolyPatch::neighbEdges() const
{
    if (!neighbEdgesPtr_.valid())
    {
        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return neighbEdgesPtr_();
}


void Foam::processorPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    if
    (
        !Pstream::parRun()
     || transform() == NOORDERING
    )
    {
        return;
    }

    if (debug)
    {
        fileName nm
        (
            boundaryMesh().mesh().time().path()
           /name()+"_faces.obj"
        );
        Pout<< "processorPolyPatch::order : Writing my " << pp.size()
            << " faces to OBJ file " << nm << endl;
        writeOBJ(nm, pp, pp.points());

        // Calculate my face centres
        const pointField& fc = pp.faceCentres();

        OFstream localStr
        (
            boundaryMesh().mesh().time().path()
           /name() + "_localFaceCentres.obj"
        );
        Pout<< "processorPolyPatch::order : "
            << "Dumping " << fc.size()
            << " local faceCentres to " << localStr.name() << endl;

        forAll(fc, facei)
        {
            writeOBJ(localStr, fc[facei]);
        }
    }

    if (owner())
    {
        if (transform() == COINCIDENTFULLMATCH)
        {
            // Pass the patch points and faces across
            UOPstream toNeighbour(neighbProcNo(), pBufs);
            toNeighbour << pp.localPoints()
                        << pp.localFaces();
        }
        else
        {
            const pointField& ppPoints = pp.points();

            pointField anchors(getAnchorPoints(pp, ppPoints, transform()));

            // Get the average of the points of each face. This is needed in
            // case the face centroid calculation is incorrect due to the face
            // having a very high aspect ratio.
            pointField facePointAverages(pp.size(), Zero);
            forAll(pp, fI)
            {
                const labelList& facePoints = pp[fI];

                forAll(facePoints, pI)
                {
                    facePointAverages[fI] += ppPoints[facePoints[pI]];
                }

                facePointAverages[fI] /= facePoints.size();
            }

            // Now send all info over to the neighbour
            UOPstream toNeighbour(neighbProcNo(), pBufs);
            toNeighbour << pp.faceCentres() << pp.faceNormals()
                        << anchors << facePointAverages;
        }
    }
}


Foam::label Foam::processorPolyPatch::matchFace
(
    const face& a,
    const pointField& aPts,
    const face& b,
    const pointField& bPts,
    const bool sameOrientation,
    const scalar absTolSqr,
    scalar& matchDistSqr
)
{
    if (a.size() != b.size())
    {
        return -1;
    }

    enum CirculatorBase::direction circulateDirection
        = CirculatorBase::CLOCKWISE;

    if (!sameOrientation)
    {
        circulateDirection = CirculatorBase::ANTICLOCKWISE;
    }

    label matchFp = -1;

    scalar closestMatchDistSqr = sqr(great);

    ConstCirculator<face> aCirc(a);
    ConstCirculator<face> bCirc(b);

    do
    {
        const scalar diffSqr = magSqr(aPts[aCirc()] - bPts[bCirc()]);

        if (diffSqr < absTolSqr)
        {
            // Found a matching point. Set the fulcrum of b to the iterator
            ConstCirculator<face> bCirc2 = bCirc;
            ++aCirc;

            bCirc2.setFulcrumToIterator();

            if (!sameOrientation)
            {
                --bCirc2;
            }
            else
            {
                ++bCirc2;
            }

            matchDistSqr = diffSqr;

            do
            {
                const scalar diffSqr2 = magSqr(aPts[aCirc()] - bPts[bCirc2()]);

                if (diffSqr2 > absTolSqr)
                {
                    // No match.
                    break;
                }

                matchDistSqr += diffSqr2;
            }
            while
            (
                aCirc.circulate(CirculatorBase::CLOCKWISE),
                bCirc2.circulate(circulateDirection)
            );

            if (!aCirc.circulate())
            {
                if (matchDistSqr < closestMatchDistSqr)
                {
                    closestMatchDistSqr = matchDistSqr;

                    if (!sameOrientation)
                    {
                        matchFp = a.size() - bCirc.nRotations();
                    }
                    else
                    {
                        matchFp = bCirc.nRotations();
                    }

                    if (closestMatchDistSqr == 0)
                    {
                        break;
                    }
                }
            }

            // Reset aCirc
            aCirc.setIteratorToFulcrum();
        }

    } while (bCirc.circulate(circulateDirection));

    matchDistSqr = closestMatchDistSqr;

    return matchFp;
}


bool Foam::processorPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Note: we only get the faces that originate from internal faces.

    if
    (
        !Pstream::parRun()
     || transform() == NOORDERING
    )
    {
        return false;
    }

    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    bool change = false;

    if (owner())
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFacei)
        {
            faceMap[patchFacei] = patchFacei;
        }

        if (transform() != COINCIDENTFULLMATCH)
        {
            const pointField& ppPoints = pp.points();

            pointField anchors(getAnchorPoints(pp, ppPoints, transform()));

            // Calculate typical distance from face centre
            scalarField tols
            (
                matchTolerance()*calcFaceTol(pp, pp.points(), pp.faceCentres())
            );

            forAll(faceMap, patchFacei)
            {
                const point& wantedAnchor = anchors[patchFacei];

                rotation[patchFacei] = getRotation
                (
                    ppPoints,
                    pp[patchFacei],
                    wantedAnchor,
                    tols[patchFacei]
                );

                if (rotation[patchFacei] > 0)
                {
                    change = true;
                }
            }
        }

        return change;
    }
    else
    {
        // Calculate the absolute matching tolerance
        scalarField tols
        (
            matchTolerance()*calcFaceTol(pp, pp.points(), pp.faceCentres())
        );

        if (transform() == COINCIDENTFULLMATCH)
        {
            vectorField masterPts;
            faceList masterFaces;

            {
                UIPstream fromNeighbour(neighbProcNo(), pBufs);
                fromNeighbour >> masterPts >> masterFaces;
            }

            const pointField& localPts = pp.localPoints();
            const faceList& localFaces = pp.localFaces();

            label nMatches = 0;

            forAll(pp, lFacei)
            {
                const face& localFace = localFaces[lFacei];
                label faceRotation = -1;

                const scalar absTolSqr = sqr(tols[lFacei]);

                scalar closestMatchDistSqr = sqr(great);
                scalar matchDistSqr = sqr(great);
                label closestFaceMatch = -1;
                label closestFaceRotation = -1;

                forAll(masterFaces, mFacei)
                {
                    const face& masterFace = masterFaces[mFacei];

                    faceRotation = matchFace
                    (
                        localFace,
                        localPts,
                        masterFace,
                        masterPts,
                        false,
                        absTolSqr,
                        matchDistSqr
                    );

                    if
                    (
                        faceRotation != -1
                     && matchDistSqr < closestMatchDistSqr
                    )
                    {
                        closestMatchDistSqr = matchDistSqr;
                        closestFaceMatch = mFacei;
                        closestFaceRotation = faceRotation;
                    }

                    if (closestMatchDistSqr == 0)
                    {
                        break;
                    }
                }

                if
                (
                    closestFaceRotation != -1
                 && closestMatchDistSqr < absTolSqr
                )
                {
                    faceMap[lFacei] = closestFaceMatch;

                    rotation[lFacei] = closestFaceRotation;

                    if (lFacei != closestFaceMatch || closestFaceRotation > 0)
                    {
                        change = true;
                    }

                    nMatches++;
                }
                else
                {
                    Pout<< "Number of matches = " << nMatches << " / "
                        << pp.size() << endl;

                    pointField pts(localFace.size());
                    forAll(localFace, pI)
                    {
                        const label localPtI = localFace[pI];
                        pts[pI] = localPts[localPtI];
                    }

                    FatalErrorInFunction
                        << "No match for face " << localFace << nl << pts
                        << abort(FatalError);
                }
            }

            return change;
        }
        else
        {
            vectorField masterCtrs;
            vectorField masterNormals;
            vectorField masterAnchors;
            vectorField masterFacePointAverages;

            // Receive data from neighbour
            {
                UIPstream fromNeighbour(neighbProcNo(), pBufs);
                fromNeighbour >> masterCtrs >> masterNormals
                              >> masterAnchors >> masterFacePointAverages;
            }

            if (debug || masterCtrs.size() != pp.size())
            {
                {
                    OFstream nbrStr
                    (
                        boundaryMesh().mesh().time().path()
                       /name() + "_nbrFaceCentres.obj"
                    );
                    Pout<< "processorPolyPatch::order : "
                        << "Dumping neighbour faceCentres to " << nbrStr.name()
                        << endl;
                    forAll(masterCtrs, facei)
                    {
                        writeOBJ(nbrStr, masterCtrs[facei]);
                    }
                }

                if (masterCtrs.size() != pp.size())
                {
                    FatalErrorInFunction
                        << "in patch:" << name() << " : "
                        << "Local size of patch is " << pp.size() << " (faces)."
                        << endl
                        << "Received from neighbour " << masterCtrs.size()
                        << " faceCentres!"
                        << abort(FatalError);
                }
            }

            // Geometric match of face centre vectors
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // 1. Try existing ordering and transformation
            bool matchedAll = matchPoints
            (
                pp.faceCentres(),
                masterCtrs,
                pp.faceNormals(),
                masterNormals,
                tols,
                false,
                faceMap
            );

            // Try using face point average for matching
            if (!matchedAll)
            {
                const pointField& ppPoints = pp.points();

                pointField facePointAverages(pp.size(), Zero);
                forAll(pp, fI)
                {
                    const labelList& facePoints = pp[fI];

                    forAll(facePoints, pI)
                    {
                        facePointAverages[fI] += ppPoints[facePoints[pI]];
                    }

                    facePointAverages[fI] /= facePoints.size();
                }

                scalarField tols2
                (
                    calcFaceTol(pp, pp.points(), facePointAverages)
                );

                labelList faceMap2(faceMap.size(), -1);
                matchedAll = matchPoints
                (
                    facePointAverages,
                    masterFacePointAverages,
                    pp.faceNormals(),
                    masterNormals,
                    tols2,
                    true,
                    faceMap2
                );

                forAll(faceMap, oldFacei)
                {
                    if (faceMap[oldFacei] == -1)
                    {
                        faceMap[oldFacei] = faceMap2[oldFacei];
                    }
                }

                matchedAll = true;

                forAll(faceMap, oldFacei)
                {
                    if (faceMap[oldFacei] == -1)
                    {
                        matchedAll = false;
                    }
                }
            }

            if (!matchedAll || debug)
            {
                // Dump faces
                fileName str
                (
                    boundaryMesh().mesh().time().path()
                   /name() + "_faces.obj"
                );
                Pout<< "processorPolyPatch::order :"
                    << " Writing faces to OBJ file " << str.name() << endl;
                writeOBJ(str, pp, pp.points());

                OFstream ccStr
                (
                    boundaryMesh().mesh().time().path()
                   /name() + "_faceCentresConnections.obj"
                );

                Pout<< "processorPolyPatch::order :"
                    << " Dumping newly found match as lines between"
                    << " corresponding face centres to OBJ file "
                    << ccStr.name()
                    << endl;

                label vertI = 0;

                forAll(pp.faceCentres(), facei)
                {
                    label masterFacei = faceMap[facei];

                    if (masterFacei != -1)
                    {
                        const point& c0 = masterCtrs[masterFacei];
                        const point& c1 = pp.faceCentres()[facei];
                        writeOBJ(ccStr, c0, c1, vertI);
                    }
                }
            }

            if (!matchedAll)
            {
                SeriousErrorInFunction
                    << "in patch:" << name() << " : "
                    << "Cannot match vectors to faces on both sides of patch"
                    << endl
                    << "    masterCtrs[0]:" << masterCtrs[0] << endl
                    << "    ctrs[0]:" << pp.faceCentres()[0] << endl
                    << "    Check your topology changes or maybe you have"
                    << " multiple separated (from cyclics) processor patches"
                    << endl
                    << "    Continuing with incorrect face ordering from now on"
                    << endl;

                return false;
            }

            // Set rotation.
            forAll(faceMap, oldFacei)
            {
                // The face f will be at newFacei (after morphing) and we want
                // its anchorPoint (= f[0]) to align with the anchorpoint for
                // the corresponding face on the other side.

                label newFacei = faceMap[oldFacei];

                const point& wantedAnchor = masterAnchors[newFacei];

                rotation[newFacei] = getRotation
                (
                    pp.points(),
                    pp[oldFacei],
                    wantedAnchor,
                    tols[oldFacei]
                );

                if (rotation[newFacei] == -1)
                {
                    SeriousErrorInFunction
                        << "in patch " << name()
                        << " : "
                        << "Cannot find point on face " << pp[oldFacei]
                        << " with vertices "
                        << UIndirectList<point>(pp.points(), pp[oldFacei])()
                        << " that matches point " << wantedAnchor
                        << " when matching the halves of processor patch "
                        << name()
                        << "Continuing with incorrect face ordering from now on"
                        << endl;

                    return false;
                }
            }

            forAll(faceMap, facei)
            {
                if (faceMap[facei] != facei || rotation[facei] != 0)
                {
                    return true;
                }
            }

            return false;
        }
    }
}


void Foam::processorPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    os.writeKeyword("myProcNo") << myProcNo_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbProcNo") << neighbProcNo_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
