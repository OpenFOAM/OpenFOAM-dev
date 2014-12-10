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

#include "pointZoneSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointZoneSet, 0);

addToRunTimeSelectionTable(topoSet, pointZoneSet, word);
addToRunTimeSelectionTable(topoSet, pointZoneSet, size);
addToRunTimeSelectionTable(topoSet, pointZoneSet, set);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pointZoneSet::updateSet()
{
    labelList order;
    sortedOrder(addressing_, order);
    inplaceReorder(order, addressing_);

    pointSet::clearStorage();
    pointSet::resize(2*addressing_.size());
    forAll(addressing_, i)
    {
        pointSet::insert(addressing_[i]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pointZoneSet::pointZoneSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    pointSet(mesh, name, 1000),  // do not read pointSet
    mesh_(mesh),
    addressing_(0)
{
    const pointZoneMesh& pointZones = mesh.pointZones();
    label zoneID = pointZones.findZoneID(name);

    if
    (
        r == IOobject::MUST_READ
     || r == IOobject::MUST_READ_IF_MODIFIED
     || (r == IOobject::READ_IF_PRESENT && zoneID != -1)
    )
    {
        const pointZone& fz = pointZones[zoneID];
        addressing_ = fz;
    }

    updateSet();

    check(mesh.nPoints());
}


pointZoneSet::pointZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    pointSet(mesh, name, size, w),
    mesh_(mesh),
    addressing_(0)
{
    updateSet();
}


pointZoneSet::pointZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    pointSet(mesh, name, set.size(), w),
    mesh_(mesh),
    addressing_(refCast<const pointZoneSet>(set).addressing())
{
    updateSet();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pointZoneSet::~pointZoneSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label pointI = 0; pointI < maxLen; pointI++)
    {
        if (!found(pointI))
        {
            n++;
        }
    }

    // Fill
    addressing_.setSize(n);
    n = 0;

    for (label pointI = 0; pointI < maxLen; pointI++)
    {
        if (!found(pointI))
        {
            addressing_[n] = pointI;
            n++;
        }
    }
    updateSet();
}


void pointZoneSet::subset(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    const pointZoneSet& fSet = refCast<const pointZoneSet>(set);

    forAll(fSet.addressing(), i)
    {
        label pointI = fSet.addressing()[i];

        if (found(pointI))
        {
            newAddressing.append(pointI);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void pointZoneSet::addSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_);

    const pointZoneSet& fSet = refCast<const pointZoneSet>(set);

    forAll(fSet.addressing(), i)
    {
        label pointI = fSet.addressing()[i];

        if (!found(pointI))
        {
            newAddressing.append(pointI);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void pointZoneSet::deleteSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    const pointZoneSet& fSet = refCast<const pointZoneSet>(set);

    forAll(addressing_, i)
    {
        label pointI = addressing_[i];

        if (!fSet.found(pointI))
        {
            // Not found in fSet so add
            newAddressing.append(pointI);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void pointZoneSet::sync(const polyMesh& mesh)
{}


label pointZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nPoints();
}


//- Write using given format, version and compression
bool pointZoneSet::writeObject
(
    IOstream::streamFormat s,
    IOstream::versionNumber v,
    IOstream::compressionType c
) const
{
    // Write shadow pointSet
    word oldTypeName = typeName;
    const_cast<word&>(type()) = pointSet::typeName;
    bool ok = pointSet::writeObject(s, v, c);
    const_cast<word&>(type()) = oldTypeName;

    // Modify pointZone
    pointZoneMesh& pointZones = const_cast<polyMesh&>(mesh_).pointZones();
    label zoneID = pointZones.findZoneID(name());

    if (zoneID == -1)
    {
        zoneID = pointZones.size();

        pointZones.setSize(zoneID+1);
        pointZones.set
        (
            zoneID,
            new pointZone
            (
                name(),
                addressing_,
                zoneID,
                pointZones
            )
        );
    }
    else
    {
        pointZones[zoneID] = addressing_;
    }
    pointZones.clearAddressing();

    return ok && pointZones.write();
}


void pointZoneSet::updateMesh(const mapPolyMesh& morphMap)
{
    // pointZone
    labelList newAddressing(addressing_.size());

    label n = 0;
    forAll(addressing_, i)
    {
        label pointI = addressing_[i];
        label newPointI = morphMap.reversePointMap()[pointI];
        if (newPointI >= 0)
        {
            newAddressing[n] = newPointI;
            n++;
        }
    }
    newAddressing.setSize(n);

    addressing_.transfer(newAddressing);

    updateSet();
}


void pointZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    pointSet::writeDebug(os, mesh, maxLen);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
