/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "faceZoneSet.H"
#include "polyTopoChangeMap.H"
#include "polyMesh.H"
#include "setToFaceZone.H"
#include "setsToFaceZone.H"
#include "syncTools.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faceZoneSet, 0);

addToRunTimeSelectionTable(topoSet, faceZoneSet, word);
addToRunTimeSelectionTable(topoSet, faceZoneSet, size);
addToRunTimeSelectionTable(topoSet, faceZoneSet, set);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void faceZoneSet::updateSet()
{
    labelList order;
    sortedOrder(addressing_, order);
    addressing_ = UIndirectList<label>(addressing_, order)();
    flipMap_ = UIndirectList<bool>(flipMap_, order)();

    faceSet::clearStorage();
    faceSet::resize(2*addressing_.size());
    forAll(addressing_, i)
    {
        faceSet::insert(addressing_[i]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

faceZoneSet::faceZoneSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    faceSet(mesh, name, 1000),  // do not read faceSet
    mesh_(mesh),
    addressing_(0),
    flipMap_(0)
{
    const faceZones& faceZones = mesh.faceZones();
    label zoneID = faceZones.findIndex(name);

    if
    (
        (r == IOobject::MUST_READ)
     || (r == IOobject::MUST_READ_IF_MODIFIED)
     || (r == IOobject::READ_IF_PRESENT && zoneID != -1)
    )
    {
        const faceZone& fz = faceZones[zoneID];
        addressing_ = fz;
        flipMap_ = fz.flipMap();
    }

    updateSet();

    check(mesh.nFaces());
}


faceZoneSet::faceZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    faceSet(mesh, name, size, w),
    mesh_(mesh),
    addressing_(0),
    flipMap_(0)
{
    updateSet();
}


faceZoneSet::faceZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    faceSet(mesh, name, set.size(), w),
    mesh_(mesh),
    addressing_(refCast<const faceZoneSet>(set).addressing()),
    flipMap_(refCast<const faceZoneSet>(set).flipMap())
{
    updateSet();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faceZoneSet::~faceZoneSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void faceZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label facei = 0; facei < maxLen; facei++)
    {
        if (!found(facei))
        {
            n++;
        }
    }

    // Fill
    addressing_.setSize(n);
    flipMap_.setSize(n);
    n = 0;

    for (label facei = 0; facei < maxLen; facei++)
    {
        if (!found(facei))
        {
            addressing_[n] = facei;
            flipMap_[n] = false;         //? or true?
            n++;
        }
    }
    updateSet();
}


void faceZoneSet::subset(const topoSet& set)
{
    label nConflict = 0;

    DynamicList<label> newAddressing(addressing_.size());
    DynamicList<bool> newFlipMap(flipMap_.size());

    Map<label> faceToIndex(addressing_.size());
    forAll(addressing_, i)
    {
        faceToIndex.insert(addressing_[i], i);
    }

    const faceZoneSet& fSet = refCast<const faceZoneSet>(set);

    forAll(fSet.addressing(), i)
    {
        label facei = fSet.addressing()[i];

        Map<label>::const_iterator iter = faceToIndex.find(facei);

        if (iter != faceToIndex.end())
        {
            label index = iter();

            if (fSet.flipMap()[i] != flipMap_[index])
            {
                nConflict++;
            }
            newAddressing.append(facei);
            newFlipMap.append(flipMap_[index]);
        }
    }

    if (nConflict > 0)
    {
        WarningInFunction
            << "subset : there are " << nConflict
            << " faces with different orientation in faceZonesSets "
            << name() << " and " << set.name() << endl;
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


void faceZoneSet::addSet(const topoSet& set)
{
    label nConflict = 0;

    DynamicList<label> newAddressing(addressing_);
    DynamicList<bool> newFlipMap(flipMap_);

    Map<label> faceToIndex(addressing_.size());
    forAll(addressing_, i)
    {
        faceToIndex.insert(addressing_[i], i);
    }

    const faceZoneSet& fSet = refCast<const faceZoneSet>(set);

    forAll(fSet.addressing(), i)
    {
        label facei = fSet.addressing()[i];

        Map<label>::const_iterator iter = faceToIndex.find(facei);

        if (iter != faceToIndex.end())
        {
            label index = iter();

            if (fSet.flipMap()[i] != flipMap_[index])
            {
                nConflict++;
            }
        }
        else
        {
            newAddressing.append(facei);
            newFlipMap.append(fSet.flipMap()[i]);
        }
    }

    if (nConflict > 0)
    {
        WarningInFunction
            << "addSet : there are " << nConflict
            << " faces with different orientation in faceZonesSets "
            << name() << " and " << set.name() << endl;
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


void faceZoneSet::deleteSet(const topoSet& set)
{
    label nConflict = 0;

    DynamicList<label> newAddressing(addressing_.size());
    DynamicList<bool> newFlipMap(flipMap_.size());

    const faceZoneSet& fSet = refCast<const faceZoneSet>(set);

    Map<label> faceToIndex(fSet.addressing().size());
    forAll(fSet.addressing(), i)
    {
        faceToIndex.insert(fSet.addressing()[i], i);
    }

    forAll(addressing_, i)
    {
        label facei = addressing_[i];

        Map<label>::const_iterator iter = faceToIndex.find(facei);

        if (iter != faceToIndex.end())
        {
            label index = iter();

            if (fSet.flipMap()[index] != flipMap_[i])
            {
                nConflict++;
            }
        }
        else
        {
            // Not found in fSet so add
            newAddressing.append(facei);
            newFlipMap.append(fSet.flipMap()[i]);
        }
    }

    if (nConflict > 0)
    {
        WarningInFunction
            << "deleteSet : there are " << nConflict
            << " faces with different orientation in faceZonesSets "
            << name() << " and " << set.name() << endl;
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


void faceZoneSet::sync(const polyMesh& mesh)
{
    // Make sure that the faceZone is consistent with the faceSet
    {
        const labelHashSet zoneSet(addressing_);

        // Get elements that are in zone but not faceSet
        labelHashSet badSet(zoneSet);
        badSet -= *this;

        // Add elements that are in faceSet but not in zone
        labelHashSet fSet(*this);
        fSet -= zoneSet;

        badSet += fSet;

        label nBad = returnReduce(badSet.size(), sumOp<label>());

        if (nBad)
        {
            WarningInFunction << "Detected " << nBad
                << " faces that are in the faceZone but not"
                << " in the faceSet or vice versa."
                << " The faceZoneSet should only be manipulated"
                << " using " << setsToFaceZone::typeName
                << " or " << setToFaceZone::typeName << endl;
        }
    }


    // Make sure that on coupled faces orientation is opposite. Pushes
    // master orientation to slave in case of conflict.


    // 0 : not in faceZone
    // 1 : in faceZone and unflipped
    //-1 : in faceZone and flipped
    const label UNFLIPPED = 1;
    const label FLIPPED = -1;
    labelList myZoneFace(mesh.nFaces()-mesh.nInternalFaces(), 0);

    forAll(addressing_, i)
    {
        label bFacei = addressing_[i]-mesh.nInternalFaces();

        if (bFacei >= 0)
        {
            if (flipMap_[i])
            {
                myZoneFace[bFacei] = FLIPPED;
            }
            else
            {
                myZoneFace[bFacei] = UNFLIPPED;
            }
        }
    }

    labelList neiZoneFace(myZoneFace);
    syncTools::swapBoundaryFaceList(mesh, neiZoneFace);


    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));


    // Rebuild faceZone addressing and flipMap
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<label> newAddressing(addressing_.size());
    DynamicList<bool> newFlipMap(flipMap_.size());

    forAll(addressing_, i)
    {
        label facei = addressing_[i];
        if (facei < mesh.nInternalFaces())
        {
            newAddressing.append(facei);
            newFlipMap.append(flipMap_[i]);
        }
    }

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        label myStat = myZoneFace[facei-mesh.nInternalFaces()];
        label neiStat = neiZoneFace[facei-mesh.nInternalFaces()];

        if (myStat == 0)
        {
            if (neiStat == UNFLIPPED)
            {
                // Neighbour is unflipped so I am flipped
                newAddressing.append(facei);
                newFlipMap.append(true);
            }
            else if (neiStat == FLIPPED)
            {
                newAddressing.append(facei);
                newFlipMap.append(false);
            }
        }
        else
        {
            if (myStat == neiStat)
            {
                // Conflict. masterFace wins
                newAddressing.append(facei);
                if (isMasterFace[facei])
                {
                    newFlipMap.append(myStat == FLIPPED);
                }
                else
                {
                    newFlipMap.append(neiStat == UNFLIPPED);
                }
            }
            else
            {
                newAddressing.append(facei);
                newFlipMap.append(myStat == FLIPPED);
            }
        }
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


label faceZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


bool faceZoneSet::writeObject
(
    IOstream::streamFormat s,
    IOstream::versionNumber v,
    IOstream::compressionType c,
    const bool write
) const
{
    // Write shadow faceSet
    word oldTypeName = typeName;
    const_cast<word&>(type()) = faceSet::typeName;
    bool ok = faceSet::writeObject(s, v, c, write);
    const_cast<word&>(type()) = oldTypeName;

    // Modify faceZone
    faceZones& faceZones = const_cast<polyMesh&>(mesh_).faceZones();
    label zoneID = faceZones.findIndex(name());

    if (zoneID == -1)
    {
        zoneID = faceZones.size();

        faceZones.setSize(zoneID+1);
        faceZones.set
        (
            zoneID,
            new faceZone
            (
                name(),
                addressing_,
                flipMap_,
                faceZones
            )
        );
    }
    else
    {
        faceZones[zoneID].resetAddressing(addressing_, flipMap_);
    }

    return ok && faceZones.write(write);
}


void faceZoneSet::topoChange(const polyTopoChangeMap& map)
{
    // faceZone
    labelList newAddressing(addressing_.size());
    boolList newFlipMap(flipMap_.size());

    label n = 0;
    forAll(addressing_, i)
    {
        label facei = addressing_[i];
        label newFacei = map.reverseFaceMap()[facei];
        if (newFacei >= 0)
        {
            newAddressing[n] = newFacei;
            newFlipMap[n] = flipMap_[i];
            n++;
        }
    }
    newAddressing.setSize(n);
    newFlipMap.setSize(n);

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);

    updateSet();
}


void faceZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    faceSet::writeDebug(os, mesh, maxLen);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
