/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    const word& patchType
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
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
        patchType
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
    myProcNo_(dict.lookup<label>("myProcNo")),
    neighbProcNo_(dict.lookup<label>("neighbProcNo")),
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
    nbrPointsPtr_.clear();
    nbrEdgesPtr_.clear();
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


void Foam::processorPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
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
            scalar magSf = magFaceAreas()[facei];
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
                const fileName patchOBJName
                (
                    boundaryMesh().mesh().time().path()/name() + "_faces.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry : Writing my "
                    << size() << " faces to " << patchOBJName << endl;

                writeOBJ(patchOBJName, *this);

                const fileName centresOBJName
                (
                    boundaryMesh().mesh().time().path()/name()
                  + "_faceCentresConnections.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry :"
                    << " Dumping lines between corresponding face centres to "
                    << centresOBJName.name() << endl;

                writeOBJ(centresOBJName, neighbFaceCentres_, faceCentres());

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
    }
}


void Foam::processorPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    processorPolyPatch::initCalcGeometry(pBufs);
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

    nbrPointsPtr_.clear();
    nbrEdgesPtr_.clear();

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

        nbrPointsPtr_.reset(new labelList(nPoints(), -1));
        labelList& nbrPoints = nbrPointsPtr_();

        forAll(nbrPointFace, nbrPointi)
        {
            // Find face and index in face on this side.
            const face& f = localFaces()[nbrPointFace[nbrPointi]];

            label index = (f.size() - nbrPointIndex[nbrPointi]) % f.size();
            label patchPointi = f[index];

            if (nbrPoints[patchPointi] == -1)
            {
                // First reference of point
                nbrPoints[patchPointi] = nbrPointi;
            }
            else if (nbrPoints[patchPointi] >= 0)
            {
                // Point already visited. Mark as duplicate.
                nbrPoints[patchPointi] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(nbrPoints, patchPointi)
        {
            if (nbrPoints[patchPointi] == -2)
            {
                nbrPoints[patchPointi] = -1;
            }
        }

        // Convert edges.
        // ~~~~~~~~~~~~~~

        nbrEdgesPtr_.reset(new labelList(nEdges(), -1));
        labelList& nbrEdges = nbrEdgesPtr_();

        forAll(nbrEdgeFace, nbrEdgeI)
        {
            // Find face and index in face on this side.
            const labelList& f = faceEdges()[nbrEdgeFace[nbrEdgeI]];
            label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
            label patchEdgeI = f[index];

            if (nbrEdges[patchEdgeI] == -1)
            {
                // First reference of edge
                nbrEdges[patchEdgeI] = nbrEdgeI;
            }
            else if (nbrEdges[patchEdgeI] >= 0)
            {
                // Edge already visited. Mark as duplicate.
                nbrEdges[patchEdgeI] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(nbrEdges, patchEdgeI)
        {
            if (nbrEdges[patchEdgeI] == -2)
            {
                nbrEdges[patchEdgeI] = -1;
            }
        }

        // Remove any addressing used for shared points/edges calculation
        // since mostly not needed.
        primitivePatch::clearOut();
    }
}


const Foam::labelList& Foam::processorPolyPatch::nbrPoints() const
{
    if (!nbrPointsPtr_.valid())
    {
        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return nbrPointsPtr_();
}


const Foam::labelList& Foam::processorPolyPatch::nbrEdges() const
{
    if (!nbrEdgesPtr_.valid())
    {
        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return nbrEdgesPtr_();
}


void Foam::processorPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (owner())
    {
        ownToNbrOrderData ownToNbr;
        autoPtr<ownToNbrDebugOrderData> ownToNbrDebugPtr
        (
            coupledPolyPatch::debug
          ? new ownToNbrDebugOrderData()
          : nullptr
        );

        coupledPolyPatch::initOrder
        (
            ownToNbr,
            ownToNbrDebugPtr,
            pp
        );

        UOPstream toNeighbour(neighbProcNo(), pBufs);
        toNeighbour << ownToNbr;
        if (coupledPolyPatch::debug)
        {
            toNeighbour << ownToNbrDebugPtr();
        }
    }
}


bool Foam::processorPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    ownToNbrOrderData ownToNbr;
    autoPtr<ownToNbrDebugOrderData> ownToNbrDebugPtr
    (
        coupledPolyPatch::debug
      ? new ownToNbrDebugOrderData()
      : nullptr
    );

    if (!owner())
    {
        UIPstream fromOwner(neighbProcNo(), pBufs);
        fromOwner >> ownToNbr;
        ownToNbr.transform(transform());
        if (coupledPolyPatch::debug)
        {
            fromOwner >> ownToNbrDebugPtr();
            ownToNbrDebugPtr->transform(transform());
        }
    }

    return
        coupledPolyPatch::order
        (
            ownToNbr,
            ownToNbrDebugPtr,
            pp,
            faceMap,
            rotation
        );
}


void Foam::processorPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    writeEntry(os, "myProcNo", myProcNo_);
    writeEntry(os, "neighbProcNo", neighbProcNo_);
}


// ************************************************************************* //
