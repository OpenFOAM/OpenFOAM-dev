/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "processorCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorCyclicPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorCyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const word& referPatchName,
    const transformType transform,
    const word& patchType
)
:
    processorPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        myProcNo,
        neighbProcNo,
        transform,
        patchType
    ),
    referPatchName_(referPatchName),
    tag_(-1),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    processorPolyPatch(name, dict, index, bm, patchType),
    referPatchName_(dict.lookup("referPatch")),
    tag_(dict.lookupOrDefault<int>("tag", -1)),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    processorPolyPatch(pp, bm),
    referPatchName_(pp.referPatchName()),
    tag_(pp.tag()),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    processorPolyPatch(pp, bm, index, newSize, newStart),
    referPatchName_(pp.referPatchName_),
    tag_(pp.tag()),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& referPatchName
)
:
    processorPolyPatch(pp, bm, index, newSize, newStart),
    referPatchName_(referPatchName),
    tag_(-1),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    processorPolyPatch(pp, bm, index, mapAddressing, newStart),
    referPatchName_(pp.referPatchName()),
    tag_(-1),
    referPatchID_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorCyclicPolyPatch::~processorCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::processorCyclicPolyPatch::tag() const
{
    if (tag_ == -1)
    {
        // Get unique tag to use for all comms. Make sure that both sides
        // use the same tag
        const cyclicPolyPatch& cycPatch = refCast<const cyclicPolyPatch>
        (
            referPatch()
        );

        if (owner())
        {
            tag_ = Hash<word>()(cycPatch.name()) % 32768u;
        }
        else
        {
            tag_ = Hash<word>()(cycPatch.neighbPatch().name()) % 32768u;
        }

        if (tag_ == Pstream::msgType() || tag_ == -1)
        {
            FatalErrorIn("processorCyclicPolyPatch::tag() const")
                << "Tag calculated from cyclic patch name " << tag_
                << " is the same as the current message type "
                << Pstream::msgType() << " or -1" << nl
                << "Please set a non-conflicting, unique, tag by hand"
                << " using the 'tag' entry"
                << exit(FatalError);
        }
        if (debug)
        {
            Pout<< "processorCyclicPolyPatch " << name() << " uses tag " << tag_
                << endl;
        }
    }
    return tag_;
}


void Foam::processorCyclicPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // Send over processorPolyPatch data
    processorPolyPatch::initGeometry(pBufs);
}


void Foam::processorCyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    // Receive and initialise processorPolyPatch data
    processorPolyPatch::calcGeometry(pBufs);

    if (Pstream::parRun())
    {

        // Where do we store the calculated transformation?
        // - on the processor patch?
        // - on the underlying cyclic patch?
        // - or do we not auto-calculate the transformation but
        //   have option of reading it.

        // Update underlying cyclic halves. Need to do both since only one
        // half might be present as a processorCyclic.
        coupledPolyPatch& pp = const_cast<coupledPolyPatch&>(referPatch());
        pp.calcGeometry
        (
            *this,
            faceCentres(),
            faceAreas(),
            faceCellCentres(),
            neighbFaceCentres(),
            neighbFaceAreas(),
            neighbFaceCellCentres()
        );

        if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(pp);
            const_cast<cyclicPolyPatch&>(cpp.neighbPatch()).calcGeometry
            (
                *this,
                neighbFaceCentres(),
                neighbFaceAreas(),
                neighbFaceCellCentres(),
                faceCentres(),
                faceAreas(),
                faceCellCentres()
            );
        }
    }
}


void Foam::processorCyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    // Recalculate geometry
    initGeometry(pBufs);
}


void Foam::processorCyclicPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField&
)
{
    calcGeometry(pBufs);
}


void Foam::processorCyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    processorPolyPatch::initUpdateMesh(pBufs);
}


void Foam::processorCyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
     referPatchID_ = -1;
     processorPolyPatch::updateMesh(pBufs);
}


void Foam::processorCyclicPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    // For now use the same algorithm as processorPolyPatch
    processorPolyPatch::initOrder(pBufs, pp);
}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::processorCyclicPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    // For now use the same algorithm as processorPolyPatch
    return processorPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::processorCyclicPolyPatch::write(Ostream& os) const
{
    processorPolyPatch::write(os);
    os.writeKeyword("referPatch") << referPatchName_
        << token::END_STATEMENT << nl;
    if (tag_ != -1)
    {
        os.writeKeyword("tag") << tag_
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
