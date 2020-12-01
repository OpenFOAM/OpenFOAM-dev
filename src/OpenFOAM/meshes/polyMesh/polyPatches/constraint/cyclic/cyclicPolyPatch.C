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

#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "matchPoints.H"
#include "EdgeMap.H"
#include "Time.H"
#include "transformField.H"
#include "SubField.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::cyclicPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initCalcGeometry(pBufs);
}


void Foam::cyclicPolyPatch::initCalcGeometry
(
    const primitivePatch& referPatch,
    pointField& nbrCtrs,
    vectorField& nbrAreas,
    pointField& nbrCc
)
{}


void Foam::cyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    static_cast<cyclicTransform&>(*this) =
        cyclicTransform
        (
            name(),
            faceCentres(),
            faceAreas(),
            *this,
            nbrPatchName(),
            nbrPatch().faceCentres(),
            nbrPatch().faceAreas(),
            nbrPatch(),
            matchTolerance()
        );
}


void Foam::cyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::cyclicPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
}


void Foam::cyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


void Foam::cyclicPolyPatch::rename(const wordList& newNames)
{
    polyPatch::rename(newNames);
    nbrPatch().nbrPatchName_ = newNames[index()];
}


void Foam::cyclicPolyPatch::reorder(const labelUList& newToOldIndex)
{
    polyPatch::reorder(newToOldIndex);
    if (nbrPatchID_ != -1)
    {
        nbrPatchID_ = findIndex(newToOldIndex, nbrPatchID_);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
    cyclicTransform(false),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr),
    ownToNbrOrderDataPtr_(nullptr),
    ownToNbrCyclicOrderDataPtr_(nullptr),
    ownToNbrDebugOrderDataPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& nbrPatchName
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
    cyclicTransform(false),
    nbrPatchName_(nbrPatchName),
    nbrPatchID_(-1),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr),
    ownToNbrOrderDataPtr_(nullptr),
    ownToNbrCyclicOrderDataPtr_(nullptr),
    ownToNbrDebugOrderDataPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    cyclicTransform(dict, false),
    nbrPatchName_(dict.lookupOrDefault("neighbourPatch", word::null)),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr),
    ownToNbrOrderDataPtr_(nullptr),
    ownToNbrCyclicOrderDataPtr_(nullptr),
    ownToNbrDebugOrderDataPtr_(nullptr)
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "No \"neighbourPatch\" provided." << endl
            << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr),
    ownToNbrOrderDataPtr_(nullptr),
    ownToNbrCyclicOrderDataPtr_(nullptr),
    ownToNbrDebugOrderDataPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& neiName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    cyclicTransform(pp),
    nbrPatchName_(neiName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr),
    ownToNbrOrderDataPtr_(nullptr),
    ownToNbrCyclicOrderDataPtr_(nullptr),
    ownToNbrDebugOrderDataPtr_(nullptr)
{
    if (neiName == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << neiName
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr),
    ownToNbrOrderDataPtr_(nullptr),
    ownToNbrCyclicOrderDataPtr_(nullptr),
    ownToNbrDebugOrderDataPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::cyclicPolyPatch::nbrPatchName() const
{
    if (nbrPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        nbrPatchName_ = boundaryMesh()[patchID].name();
    }
    return nbrPatchName_;
}


Foam::label Foam::cyclicPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << nbrPatchName()
                << endl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic
        const cyclicPolyPatch& nbrPatch = refCast<const cyclicPolyPatch>
        (
            this->boundaryMesh()[nbrPatchID_]
        );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << endl << " but that in return specifies "
                << nbrPatch.nbrPatchName()
                << endl;
        }
    }
    return nbrPatchID_;
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledPoints() const
{
    if (!coupledPointsPtr_)
    {
        const faceList& nbrLocalFaces = nbrPatch().localFaces();
        const labelList& nbrMeshPoints = nbrPatch().meshPoints();

        // Now all we know is that relative face index in *this is same
        // as coupled face in nbrPatch and also that the 0th vertex
        // corresponds.

        // From local point to nbrPatch or -1.
        labelList coupledPoint(nPoints(), -1);

        forAll(*this, patchFacei)
        {
            const face& fA = localFaces()[patchFacei];
            const face& fB = nbrLocalFaces[patchFacei];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (coupledPoint[patchPointA] == -1)
                {
                    label indexB = (fB.size() - indexA) % fB.size();

                    // Filter out points on wedge axis
                    if (meshPoints()[patchPointA] != nbrMeshPoints[fB[indexB]])
                    {
                        coupledPoint[patchPointA] = fB[indexB];
                    }
                }
            }
        }

        coupledPointsPtr_ = new edgeList(nPoints());
        edgeList& connected = *coupledPointsPtr_;

        // Extract coupled points.
        label connectedI = 0;

        forAll(coupledPoint, i)
        {
            if (coupledPoint[i] != -1)
            {
                connected[connectedI++] = edge(i, coupledPoint[i]);
            }
        }

        connected.setSize(connectedI);

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /name() + "_coupledPoints.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with coordinates of "
                << "coupled points" << endl;

            forAll(connected, i)
            {
                const point& a = points()[meshPoints()[connected[i][0]]];
                const point& b = points()[nbrMeshPoints[connected[i][1]]];

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledPointsPtr_;
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledEdges() const
{
    if (!coupledEdgesPtr_)
    {
        const edgeList& pointCouples = coupledPoints();

        // Build map from points on *this (A) to points on neighbourpatch (B)
        Map<label> aToB(2*pointCouples.size());

        forAll(pointCouples, i)
        {
            const edge& e = pointCouples[i];

            aToB.insert(e[0], e[1]);
        }

        // Map from edge on A to points (in B indices)
        EdgeMap<label> edgeMap(nEdges());

        forAll(*this, patchFacei)
        {
            const labelList& fEdges = faceEdges()[patchFacei];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Convert edge end points to corresponding points on B side.
                Map<label>::const_iterator fnd0 = aToB.find(e[0]);
                if (fnd0 != aToB.end())
                {
                    Map<label>::const_iterator fnd1 = aToB.find(e[1]);
                    if (fnd1 != aToB.end())
                    {
                        edgeMap.insert(edge(fnd0(), fnd1()), edgeI);
                    }
                }
            }
        }

        // Use the edgeMap to get the edges on the B side.

        const cyclicPolyPatch& nbrPatch = this->nbrPatch();
        const labelList& nbrMp = nbrPatch.meshPoints();
        const labelList& mp = meshPoints();



        coupledEdgesPtr_ = new edgeList(edgeMap.size());
        edgeList& coupledEdges = *coupledEdgesPtr_;
        label coupleI = 0;

        forAll(nbrPatch, patchFacei)
        {
            const labelList& fEdges = nbrPatch.faceEdges()[patchFacei];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = nbrPatch.edges()[edgeI];

                // Look up A edge from HashTable.
                EdgeMap<label>::iterator iter = edgeMap.find(e);

                if (iter != edgeMap.end())
                {
                    label edgeA = iter();
                    const edge& eA = edges()[edgeA];

                    // Store correspondence. Filter out edges on wedge axis.
                    if
                    (
                        edge(mp[eA[0]], mp[eA[1]])
                     != edge(nbrMp[e[0]], nbrMp[e[1]])
                    )
                    {
                        coupledEdges[coupleI++] = edge(edgeA, edgeI);
                    }

                    // Remove so we build unique list
                    edgeMap.erase(iter);
                }
            }
        }
        coupledEdges.setSize(coupleI);


        // Some checks

        forAll(coupledEdges, i)
        {
            const edge& e = coupledEdges[i];

            if (e[0] < 0 || e[1] < 0)
            {
                FatalErrorInFunction
                    << "Problem : at position " << i
                    << " illegal couple:" << e
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /name() + "_coupledEdges.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with centres of "
                << "coupled edges" << endl;

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                const point& a = edges()[e[0]].centre(localPoints());
                const point& b = nbrPatch.edges()[e[1]].centre
                (
                    nbrPatch.localPoints()
                );

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledEdgesPtr_;
}


void Foam::cyclicPolyPatch::initOrder
(
    PstreamBuffers&,
    const primitivePatch& pp
) const
{
    if (pp.empty())
    {
        return;
    }

    if (owner())
    {
        ownToNbrOrderDataPtr_ = new ownToNbrOrderData();
        if (coupledPolyPatch::debug)
        {
            ownToNbrDebugOrderDataPtr_ = new ownToNbrDebugOrderData();
        }

        coupledPolyPatch::initOrder
        (
            ownToNbrOrderDataPtr_(),
            ownToNbrDebugOrderDataPtr_,
            pp
        );

        const scalarField magAreas(mag(pp.faceAreas()));

        ownToNbrCyclicOrderDataPtr_ = new ownToNbrCyclicOrderData();
        ownToNbrCyclicOrderDataPtr_->ctr =
            sum(pp.faceCentres()*magAreas)/sum(magAreas);
        ownToNbrCyclicOrderDataPtr_->area = sum(pp.faceAreas());
    }
}


bool Foam::cyclicPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (pp.empty())
    {
        return false;
    }

    ownToNbrOrderData ownToNbr;
    autoPtr<ownToNbrDebugOrderData> ownToNbrDebugPtr(nullptr);

    if (!owner())
    {
        ownToNbr = nbrPatch().ownToNbrOrderDataPtr_();
        ownToNbrDebugPtr = nbrPatch().ownToNbrDebugOrderDataPtr_;

        cyclicTransform ct
        (
            name(),
            pp.faceCentres(),
            pp.faceAreas(),
            *this,
            nbrPatchName(),
            pointField(1, nbrPatch().ownToNbrCyclicOrderDataPtr_->ctr),
            vectorField(1, nbrPatch().ownToNbrCyclicOrderDataPtr_->area),
            nbrPatch(),
            matchTolerance()
        );

        ownToNbr.transform(ct.transform());
        if (ownToNbrDebugPtr.valid())
        {
            ownToNbrDebugPtr->transform(ct.transform());
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


void Foam::cyclicPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);

    if (!nbrPatchName_.empty())
    {
        writeEntry(os, "neighbourPatch", nbrPatchName_);
    }

    coupleGroup_.write(os);

    cyclicTransform::write(os);
}


// ************************************************************************* //
