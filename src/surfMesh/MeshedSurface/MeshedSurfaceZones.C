/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "MeshedSurface.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::checkZones()
{
    // extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently
    surfZoneList& zones = this->storedZones();
    if (zones.size())
    {
        label count = 0;
        forAll(zones, zoneI)
        {
            zones[zoneI].start() = count;
            count += zones[zoneI].size();
        }

        if (count < this->size())
        {
            WarningInFunction
                << "more faces " << this->size() << " than zones " << count
                << " ... extending final zone"
                << endl;

            zones.last().size() += count - this->size();
        }
        else if (count > this->size())
        {
            FatalErrorInFunction
                << "more zones " << count << " than faces " << this->size()
                << exit(FatalError);
        }
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::sortFacesAndStore
(
    List<Face>&& unsortedFaces,
    List<label>&& zoneIds,
    const bool sorted
)
{
    List<Face>  oldFaces(unsortedFaces);
    List<label> zones(zoneIds);

    if (sorted)
    {
        // already sorted - simply transfer faces
        this->storedFaces().transfer(oldFaces);
    }
    else
    {
        // unsorted - determine the sorted order:
        // avoid SortableList since we discard the main list anyhow
        List<label> faceMap;
        sortedOrder(zones, faceMap);
        zones.clear();

        // sorted faces
        List<Face> newFaces(faceMap.size());
        forAll(faceMap, facei)
        {
            // use transfer to recover memory where possible
            newFaces[facei].transfer(oldFaces[faceMap[facei]]);
        }
        this->storedFaces().transfer(newFaces);
    }
    zones.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const UList<surfZone>& srfZones,
    const bool cullEmpty
)
{
    label nZone = 0;

    surfZoneList& zones = this->storedZones();
    zones.setSize(zones.size());
    forAll(zones, zoneI)
    {
        if (srfZones[zoneI].size() || !cullEmpty)
        {
            zones[nZone] = surfZone(srfZones[zoneI], nZone);
            nZone++;
        }
    }
    zones.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const labelUList& sizes,
    const UList<word>& names,
    const bool cullEmpty
)
{
    label start   = 0;
    label nZone = 0;

    surfZoneList& zones = this->storedZones();
    zones.setSize(sizes.size());
    forAll(zones, zoneI)
    {
        if (sizes[zoneI] || !cullEmpty)
        {
            zones[nZone] = surfZone
            (
                names[zoneI],
                sizes[zoneI],
                start,
                nZone
            );
            start += sizes[zoneI];
            nZone++;
        }
    }
    zones.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const labelUList& sizes,
    const bool cullEmpty
)
{
    label start   = 0;
    label nZone = 0;

    surfZoneList& zones = this->storedZones();
    zones.setSize(sizes.size());
    forAll(zones, zoneI)
    {
        if (sizes[zoneI] || !cullEmpty)
        {
            zones[nZone] = surfZone
            (
                word("zone") + ::Foam::name(nZone),
                sizes[zoneI],
                start,
                nZone
            );
            start += sizes[zoneI];
            nZone++;
        }
    }
    zones.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::removeZones()
{
    this->storedZones().clear();
}


// ************************************************************************* //
