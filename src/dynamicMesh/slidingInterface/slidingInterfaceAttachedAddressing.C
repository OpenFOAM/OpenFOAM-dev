/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "slidingInterface.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::calcAttachedAddressing() const
{
    if (debug)
    {
        Pout<< "void Foam::slidingInterface::calcAttachedAddressing() const "
            << " for object " << name() << " : "
            << "Calculating zone face-cell addressing."
            << endl;
    }

    if (!attached_)
    {
        // Clear existing addressing
        clearAttachedAddressing();

        const polyMesh& mesh = topoChanger().mesh();
        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();
        const faceZoneMesh& faceZones = mesh.faceZones();

        // Master side

        const primitiveFacePatch& masterPatch =
            faceZones[masterFaceZoneID_.index()]();

        const labelList& masterPatchFaces =
            faceZones[masterFaceZoneID_.index()];

        const boolList& masterFlip =
            faceZones[masterFaceZoneID_.index()].flipMap();

        masterFaceCellsPtr_ = new labelList(masterPatchFaces.size());
        labelList& mfc = *masterFaceCellsPtr_;

        forAll(masterPatchFaces, faceI)
        {
            if (masterFlip[faceI])
            {
                mfc[faceI] = nei[masterPatchFaces[faceI]];
            }
            else
            {
                mfc[faceI] = own[masterPatchFaces[faceI]];
            }
        }

        // Slave side

        const primitiveFacePatch& slavePatch =
            faceZones[slaveFaceZoneID_.index()]();

        const labelList& slavePatchFaces =
            faceZones[slaveFaceZoneID_.index()];

        const boolList& slaveFlip =
            faceZones[slaveFaceZoneID_.index()].flipMap();

        slaveFaceCellsPtr_ = new labelList(slavePatchFaces.size());
        labelList& sfc = *slaveFaceCellsPtr_;

        forAll(slavePatchFaces, faceI)
        {
            if (slaveFlip[faceI])
            {
                sfc[faceI] = nei[slavePatchFaces[faceI]];
            }
            else
            {
                sfc[faceI] = own[slavePatchFaces[faceI]];
            }
        }

        // Check that the addressing is valid
        if (min(mfc) < 0 || min(sfc) < 0)
        {
            if (debug)
            {
                forAll(mfc, faceI)
                {
                    if (mfc[faceI] < 0)
                    {
                        Pout<< "No cell next to master patch face " << faceI
                            << ".  Global face no: " << mfc[faceI]
                            << " own: " << own[masterPatchFaces[faceI]]
                            << " nei: " << nei[masterPatchFaces[faceI]]
                            << " flip: " << masterFlip[faceI] << endl;
                    }
                }

                forAll(sfc, faceI)
                {
                    if (sfc[faceI] < 0)
                    {
                        Pout<< "No cell next to slave patch face " << faceI
                            << ".  Global face no: " << sfc[faceI]
                            << " own: " << own[slavePatchFaces[faceI]]
                            << " nei: " << nei[slavePatchFaces[faceI]]
                            << " flip: " << slaveFlip[faceI] << endl;
                    }
                }
            }

            FatalErrorIn
            (
                "void slidingInterface::calcAttachedAddressing()"
                "const"
            )   << "Error is zone face-cell addressing.  Probable error in "
                << "decoupled mesh or sliding interface definition."
                << abort(FatalError);
        }

        // Calculate stick-out faces
        const labelListList& pointFaces = mesh.pointFaces();

        // Master side
        labelHashSet masterStickOutFaceMap
        (
            primitiveMesh::facesPerCell_*(masterPatch.size())
        );

        const labelList& masterMeshPoints = masterPatch.meshPoints();

        forAll(masterMeshPoints, pointI)
        {
            const labelList& curFaces = pointFaces[masterMeshPoints[pointI]];

            forAll(curFaces, faceI)
            {
                // Check if the face belongs to the master face zone;
                // if not add it
                if
                (
                    faceZones.whichZone(curFaces[faceI])
                 != masterFaceZoneID_.index()
                )
                {
                    masterStickOutFaceMap.insert(curFaces[faceI]);
                }
            }
        }

        masterStickOutFacesPtr_ = new labelList(masterStickOutFaceMap.toc());

        // Slave side
        labelHashSet slaveStickOutFaceMap
        (
            primitiveMesh::facesPerCell_*(slavePatch.size())
        );

        const labelList& slaveMeshPoints = slavePatch.meshPoints();

        forAll(slaveMeshPoints, pointI)
        {
            const labelList& curFaces = pointFaces[slaveMeshPoints[pointI]];

            forAll(curFaces, faceI)
            {
                // Check if the face belongs to the slave face zone;
                // if not add it
                if
                (
                    faceZones.whichZone(curFaces[faceI])
                 != slaveFaceZoneID_.index()
                )
                {
                    slaveStickOutFaceMap.insert(curFaces[faceI]);
                }
            }
        }

        slaveStickOutFacesPtr_ = new labelList(slaveStickOutFaceMap.toc());


        // Retired point addressing does not exist at this stage.
        // It will be filled when the interface is coupled.
        retiredPointMapPtr_ =
            new Map<label>
            (
                2*faceZones[slaveFaceZoneID_.index()]().nPoints()
            );

        // Ditto for cut point edge map.  This is a rough guess of its size
        //
        cutPointEdgePairMapPtr_ =
            new Map<Pair<edge> >
            (
                faceZones[slaveFaceZoneID_.index()]().nEdges()
            );
    }
    else
    {
        FatalErrorIn
        (
            "void slidingInterface::calcAttachedAddressing() const"
        )   << "The interface is attached.  The zone face-cell addressing "
            << "cannot be assembled for object " << name()
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "void Foam::slidingInterface::calcAttachedAddressing() const "
            << " for object " << name() << " : "
            << "Finished calculating zone face-cell addressing."
            << endl;
    }
}


void Foam::slidingInterface::clearAttachedAddressing() const
{
    deleteDemandDrivenData(masterFaceCellsPtr_);
    deleteDemandDrivenData(slaveFaceCellsPtr_);

    deleteDemandDrivenData(masterStickOutFacesPtr_);
    deleteDemandDrivenData(slaveStickOutFacesPtr_);

    deleteDemandDrivenData(retiredPointMapPtr_);
    deleteDemandDrivenData(cutPointEdgePairMapPtr_);
}


void Foam::slidingInterface::renumberAttachedAddressing
(
    const mapPolyMesh& m
) const
{
    // Get reference to reverse cell renumbering
    // The renumbering map is needed the other way around, i.e. giving
    // the new cell number for every old cell next to the interface.
    const labelList& reverseCellMap = m.reverseCellMap();

    const labelList& mfc = masterFaceCells();
    const labelList& sfc = slaveFaceCells();

    // Master side
    labelList* newMfcPtr = new labelList(mfc.size(), -1);
    labelList& newMfc = *newMfcPtr;

    const labelList& mfzRenumber =
        m.faceZoneFaceMap()[masterFaceZoneID_.index()];

    forAll(mfc, faceI)
    {
        label newCellI = reverseCellMap[mfc[mfzRenumber[faceI]]];

        if (newCellI >= 0)
        {
            newMfc[faceI] = newCellI;
        }
    }

    // Slave side
    labelList* newSfcPtr = new labelList(sfc.size(), -1);
    labelList& newSfc = *newSfcPtr;

    const labelList& sfzRenumber =
        m.faceZoneFaceMap()[slaveFaceZoneID_.index()];

    forAll(sfc, faceI)
    {
        label newCellI = reverseCellMap[sfc[sfzRenumber[faceI]]];

        if (newCellI >= 0)
        {
            newSfc[faceI] = newCellI;
        }
    }

    if (debug)
    {
        // Check if all the mapped cells are live
        if (min(newMfc) < 0 || min(newSfc) < 0)
        {
            FatalErrorIn
            (
                "void slidingInterface::renumberAttachedAddressing("
                "const mapPolyMesh& m) const"
            )   << "Error in cell renumbering for object " << name()
                << ".  Some of master cells next "
                << "to the interface have been removed."
                << abort(FatalError);
        }
    }

    // Renumber stick-out faces

    const labelList& reverseFaceMap = m.reverseFaceMap();

    // Master side
    const labelList& msof = masterStickOutFaces();

    labelList* newMsofPtr = new labelList(msof.size(), -1);
    labelList& newMsof = *newMsofPtr;

    forAll(msof, faceI)
    {
        label newFaceI = reverseFaceMap[msof[faceI]];

        if (newFaceI >= 0)
        {
            newMsof[faceI] = newFaceI;
        }
    }
//     Pout<< "newMsof: " << newMsof << endl;
    // Slave side
    const labelList& ssof = slaveStickOutFaces();

    labelList* newSsofPtr = new labelList(ssof.size(), -1);
    labelList& newSsof = *newSsofPtr;

    forAll(ssof, faceI)
    {
        label newFaceI = reverseFaceMap[ssof[faceI]];

        if (newFaceI >= 0)
        {
            newSsof[faceI] = newFaceI;
        }
    }
//     Pout<< "newSsof: " << newSsof << endl;
    if (debug)
    {
        // Check if all the mapped cells are live
        if (min(newMsof) < 0 || min(newSsof) < 0)
        {
            FatalErrorIn
            (
                "void slidingInterface::renumberAttachedAddressing("
                "const mapPolyMesh& m) const"
            )   << "Error in face renumbering for object " << name()
                << ".  Some of stick-out next "
                << "to the interface have been removed."
                << abort(FatalError);
        }
    }

    // Renumber the retired point map. Need to take a copy!
    const Map<label> rpm = retiredPointMap();

    Map<label>* newRpmPtr = new Map<label>(rpm.size());
    Map<label>& newRpm = *newRpmPtr;

    const labelList rpmToc = rpm.toc();

    // Get reference to point renumbering
    const labelList& reversePointMap = m.reversePointMap();

    label key, value;

    forAll(rpmToc, rpmTocI)
    {
        key = reversePointMap[rpmToc[rpmTocI]];

        value = reversePointMap[rpm.find(rpmToc[rpmTocI])()];

        if (debug)
        {
            // Check if all the mapped cells are live
            if (key < 0 || value < 0)
            {
                FatalErrorIn
                (
                    "void slidingInterface::renumberAttachedAddressing("
                    "const mapPolyMesh& m) const"
                )   << "Error in retired point numbering for object "
                    << name() << ".  Some of master "
                    << "points have been removed."
                    << abort(FatalError);
            }
        }

        newRpm.insert(key, value);
    }

    // Renumber the cut point edge pair map. Need to take a copy!
    const Map<Pair<edge> > cpepm = cutPointEdgePairMap();

    Map<Pair<edge> >* newCpepmPtr = new Map<Pair<edge> >(cpepm.size());
    Map<Pair<edge> >& newCpepm = *newCpepmPtr;

    const labelList cpepmToc = cpepm.toc();

    forAll(cpepmToc, cpepmTocI)
    {
        key = reversePointMap[cpepmToc[cpepmTocI]];

        const Pair<edge>& oldPe = cpepm.find(cpepmToc[cpepmTocI])();

        // Re-do the edges in global addressing
        const label ms = reversePointMap[oldPe.first().start()];
        const label me = reversePointMap[oldPe.first().end()];

        const label ss = reversePointMap[oldPe.second().start()];
        const label se = reversePointMap[oldPe.second().end()];

        if (debug)
        {
            // Check if all the mapped cells are live
            if (key < 0 || ms < 0 || me < 0 || ss < 0 || se < 0)
            {
                FatalErrorIn
                (
                    "void slidingInterface::renumberAttachedAddressing("
                    "const mapPolyMesh& m) const"
                )   << "Error in cut point edge pair map numbering for object "
                    << name() << ".  Some of master points have been removed."
                    << abort(FatalError);
            }
        }

        newCpepm.insert(key, Pair<edge>(edge(ms, me), edge(ss, se)));
    }

    if (!projectedSlavePointsPtr_)
    {
        FatalErrorIn
        (
            "void slidingInterface::renumberAttachedAddressing("
            "const mapPolyMesh& m) const"
        )   << "Error in projected point numbering for object " << name()
            << abort(FatalError);
    }

    // Renumber the projected slave zone points
    const pointField& projectedSlavePoints = *projectedSlavePointsPtr_;

    pointField* newProjectedSlavePointsPtr
    (
        new pointField(projectedSlavePoints.size())
    );
    pointField& newProjectedSlavePoints = *newProjectedSlavePointsPtr;

    const labelList& sfzPointRenumber =
        m.faceZonePointMap()[slaveFaceZoneID_.index()];

    forAll(newProjectedSlavePoints, pointI)
    {
        if (sfzPointRenumber[pointI] > -1)
        {
            newProjectedSlavePoints[pointI] =
                projectedSlavePoints[sfzPointRenumber[pointI]];
        }
    }

    // Re-set the lists
    clearAttachedAddressing();

    deleteDemandDrivenData(projectedSlavePointsPtr_);

    masterFaceCellsPtr_ = newMfcPtr;
    slaveFaceCellsPtr_ = newSfcPtr;

    masterStickOutFacesPtr_ = newMsofPtr;
    slaveStickOutFacesPtr_ = newSsofPtr;

    retiredPointMapPtr_ = newRpmPtr;
    cutPointEdgePairMapPtr_ = newCpepmPtr;
    projectedSlavePointsPtr_ = newProjectedSlavePointsPtr;
}


const Foam::labelList& Foam::slidingInterface::masterFaceCells() const
{
    if (!masterFaceCellsPtr_)
    {
        FatalErrorIn
        (
            "const labelList& slidingInterface::masterFaceCells() const"
        )   << "Master zone face-cell addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *masterFaceCellsPtr_;
}


const Foam::labelList& Foam::slidingInterface::slaveFaceCells() const
{
    if (!slaveFaceCellsPtr_)
    {
        FatalErrorIn
        (
            "const labelList& slidingInterface::slaveFaceCells() const"
        )   << "Slave zone face-cell addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *slaveFaceCellsPtr_;
}


const Foam::labelList& Foam::slidingInterface::masterStickOutFaces() const
{
    if (!masterStickOutFacesPtr_)
    {
        FatalErrorIn
        (
            "const labelList& slidingInterface::masterStickOutFaces() const"
        )   << "Master zone stick-out face addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *masterStickOutFacesPtr_;
}


const Foam::labelList& Foam::slidingInterface::slaveStickOutFaces() const
{
    if (!slaveStickOutFacesPtr_)
    {
        FatalErrorIn
        (
            "const labelList& slidingInterface::slaveStickOutFaces() const"
        )   << "Slave zone stick-out face addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *slaveStickOutFacesPtr_;
}


const Foam::Map<Foam::label>& Foam::slidingInterface::retiredPointMap() const
{
    if (!retiredPointMapPtr_)
    {
        FatalErrorIn
        (
            "const Map<label>& slidingInterface::retiredPointMap() const"
        )   << "Retired point map not available for object " << name()
            << abort(FatalError);
    }

    return *retiredPointMapPtr_;
}


const Foam::Map<Foam::Pair<Foam::edge> >&
Foam::slidingInterface::cutPointEdgePairMap() const
{
    if (!cutPointEdgePairMapPtr_)
    {
        FatalErrorIn
        (
            "const Map<Pair<edge> >& slidingInterface::"
            "cutPointEdgePairMap() const"
        )   << "Retired point map not available for object " << name()
            << abort(FatalError);
    }

    return *cutPointEdgePairMapPtr_;
}


// ************************************************************************* //
