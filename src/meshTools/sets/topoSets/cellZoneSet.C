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

#include "cellZoneSet.H"
#include "polyTopoChangeMap.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cellZoneSet, 0);

addToRunTimeSelectionTable(topoSet, cellZoneSet, word);
addToRunTimeSelectionTable(topoSet, cellZoneSet, size);
addToRunTimeSelectionTable(topoSet, cellZoneSet, set);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void cellZoneSet::updateSet()
{
    labelList order;
    sortedOrder(addressing_, order);
    inplaceReorder(order, addressing_);

    cellSet::clearStorage();
    cellSet::resize(2*addressing_.size());
    forAll(addressing_, i)
    {
        cellSet::insert(addressing_[i]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    cellSet(mesh, name, 1000),  // do not read cellSet
    mesh_(mesh),
    addressing_(0)
{
    const cellZoneList& cellZones = mesh.cellZones();
    label zoneID = cellZones.findIndex(name);

    if
    (
        (r == IOobject::MUST_READ)
     || (r == IOobject::MUST_READ_IF_MODIFIED)
     || (r == IOobject::READ_IF_PRESENT && zoneID != -1)
    )
    {
        const cellZone& fz = cellZones[zoneID];
        addressing_ = fz;
    }

    updateSet();

    check(mesh.nCells());
}


cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    cellSet(mesh, name, size, w),
    mesh_(mesh),
    addressing_(0)
{
    updateSet();
}


cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    cellSet(mesh, name, set.size(), w),
    mesh_(mesh),
    addressing_(refCast<const cellZoneSet>(set).addressing())
{
    updateSet();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellZoneSet::~cellZoneSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cellZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label celli = 0; celli < maxLen; celli++)
    {
        if (!found(celli))
        {
            n++;
        }
    }

    // Fill
    addressing_.setSize(n);
    n = 0;

    for (label celli = 0; celli < maxLen; celli++)
    {
        if (!found(celli))
        {
            addressing_[n] = celli;
            n++;
        }
    }

    updateSet();
}


void cellZoneSet::subset(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    const cellZoneSet& fSet = refCast<const cellZoneSet>(set);

    forAll(fSet.addressing(), i)
    {
        label celli = fSet.addressing()[i];

        if (found(celli))
        {
            newAddressing.append(celli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void cellZoneSet::addSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_);

    const cellZoneSet& fSet = refCast<const cellZoneSet>(set);

    forAll(fSet.addressing(), i)
    {
        label celli = fSet.addressing()[i];

        if (!found(celli))
        {
            newAddressing.append(celli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void cellZoneSet::deleteSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    const cellZoneSet& fSet = refCast<const cellZoneSet>(set);

    forAll(addressing_, i)
    {
        label celli = addressing_[i];

        if (!fSet.found(celli))
        {
            // Not found in fSet so add
            newAddressing.append(celli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void cellZoneSet::sync(const polyMesh& mesh)
{
    cellSet::sync(mesh);

    // Take over contents of cellSet into addressing.
    addressing_ = sortedToc();
    updateSet();
}


label cellZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


bool cellZoneSet::writeObject
(
    IOstream::streamFormat s,
    IOstream::versionNumber v,
    IOstream::compressionType c,
    const bool write
) const
{
    // Write shadow cellSet
    word oldTypeName = typeName;
    const_cast<word&>(type()) = cellSet::typeName;
    bool ok = cellSet::writeObject(s, v, c, write);
    const_cast<word&>(type()) = oldTypeName;

    // Modify cellZone
    cellZoneList& cellZones = const_cast<polyMesh&>(mesh_).cellZones();
    label zoneID = cellZones.findIndex(name());

    if (zoneID == -1)
    {
        zoneID = cellZones.size();

        cellZones.setSize(zoneID+1);
        cellZones.set
        (
            zoneID,
            new cellZone
            (
                name(),
                addressing_,
                cellZones
            )
        );
    }
    else
    {
        cellZones[zoneID] = addressing_;
    }

    return ok && cellZones.write(write);
}


void cellZoneSet::topoChange(const polyTopoChangeMap& map)
{
    // cellZone
    labelList newAddressing(addressing_.size());

    label n = 0;
    forAll(addressing_, i)
    {
        label celli = addressing_[i];
        label newCelli = map.reverseCellMap()[celli];
        if (newCelli >= 0)
        {
            newAddressing[n] = newCelli;
            n++;
        }
    }
    newAddressing.setSize(n);

    addressing_.transfer(newAddressing);

    updateSet();
}


void cellZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    cellSet::writeDebug(os, mesh, maxLen);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
