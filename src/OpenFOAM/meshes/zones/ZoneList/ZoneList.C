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

#include "ZoneList.H"
#include "Pstream.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::read()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningInFunction
                << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        PtrListDictionary<ZoneType>& zones = *this;

        // Read zones
        Istream& is = readStream(ZonesType::typeName);

        PtrList<entry> patchEntries(is);
        zones.setSize(patchEntries.size());

        forAll(zones, zi)
        {
            zones.set
            (
                zi,
                patchEntries[zi].keyword(),
                ZoneType::New
                (
                    patchEntries[zi].keyword(),
                    patchEntries[zi].dict(),
                    static_cast<const ZonesType&>(*this)
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "ZoneList::ZoneList"
            "(const IOobject&, const MeshType&)"
        );

        close();

        return true;
    }
    else
    {
        // Nothing read
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType, class MeshType>
Foam::ZoneList<ZoneType, ZonesType, MeshType>::ZoneList
(
    const IOobject& io,
    const MeshType& mesh
)
:
    regIOobject(io),
    PtrListDictionary<ZoneType>(0),
    mesh_(mesh)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType, class MeshType>
Foam::ZoneList<ZoneType, ZonesType, MeshType>::~ZoneList()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::found
(
    const label objectIndex
) const
{
    forAll(*this, zi)
    {
        if (this->operator[](zi).localIndex(objectIndex) != -1)
        {
            return true;
        }
    }

    return false;
}


template<class ZoneType, class ZonesType, class MeshType>
Foam::labelList Foam::ZoneList<ZoneType, ZonesType, MeshType>::whichZones
(
    const label objectIndex
) const
{
    labelList zones;

    forAll(*this, zi)
    {
        if (this->operator[](zi).localIndex(objectIndex) != -1)
        {
            zones.append(zi);
        }
    }

    return zones;
}


template<class ZoneType, class ZonesType, class MeshType>
Foam::wordList Foam::ZoneList<ZoneType, ZonesType, MeshType>::types() const
{
    const PtrListDictionary<ZoneType>& zones = *this;

    wordList lst(zones.size());

    forAll(zones, zi)
    {
        lst[zi] = zones[zi].type();
    }

    return lst;
}


template<class ZoneType, class ZonesType, class MeshType>
Foam::labelHashSet Foam::ZoneList<ZoneType, ZonesType, MeshType>::zoneSet
(
    const UList<wordRe>& zoneNames,
    const bool warnNotFound
) const
{
    labelHashSet set;

    if (zoneNames.size())
    {
        forAll(zoneNames, i)
        {
            const labelList indices
            (
                this->findIndices(zoneNames[i])
            );

            if (indices.size())
            {
                set.insert(indices);
            }
            else if (warnNotFound)
            {
                WarningInFunction
                    << "Cannot find zone " << zoneNames[i]
                    << " of type " << type()
                    << endl;
            }
        }
    }

    return set;
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::append
(
    ZoneType* zonePtr
) const
{
    ZoneList<ZoneType, ZonesType, MeshType>& zones =
        const_cast<ZoneList<ZoneType, ZonesType, MeshType>&>(*this);

    if (found(zonePtr->name()))
    {
        zones[zonePtr->name()] = *zonePtr;
        delete zonePtr;
    }
    else
    {
        zones.PtrListDictionary<ZoneType>::append
        (
            zonePtr->name(),
            zonePtr
        );
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::append
(
    const ZoneType& zone
) const
{
    ZoneList<ZoneType, ZonesType, MeshType>& zones =
        const_cast<ZoneList<ZoneType, ZonesType, MeshType>&>(*this);

    if (found(zone.name()))
    {
        zones[zone.name()] = zone;
    }
    else
    {
        zones.PtrListDictionary<ZoneType>::append
        (
            zone.name(),
            zone.clone(*this)
        );
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::clearAddressing()
{
    PtrListDictionary<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].clearAddressing();
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::clear()
{
    clearAddressing();
    PtrListDictionary<ZoneType>::clear();
}


template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::checkDefinition
(
    const bool report
) const
{
    bool inError = false;

    const PtrListDictionary<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        inError |= zones[zi].checkDefinition(report);
    }
    return inError;
}


template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::checkParallelSync
(
    const bool report
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }


    const PtrListDictionary<ZoneType>& zones = *this;

    bool hasError = false;

    // Collect all names
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = this->toc();
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    List<wordList> allTypes(Pstream::nProcs());
    allTypes[Pstream::myProcNo()] = this->types();
    Pstream::gatherList(allTypes);
    Pstream::scatterList(allTypes);

    // Have every processor check but only master print error.

    for (label proci = 1; proci < allNames.size(); proci++)
    {
        if
        (
            (allNames[proci] != allNames[0])
         || (allTypes[proci] != allTypes[0])
        )
        {
            hasError = true;

            if (debug || (report && Pstream::master()))
            {
                Info<< " ***Inconsistent zones across processors, "
                       "processor 0 has zone names:" << allNames[0]
                    << " zone types:" << allTypes[0]
                    << " processor " << proci << " has zone names:"
                    << allNames[proci]
                    << " zone types:" << allTypes[proci]
                    << endl;
            }
        }
    }

    // Check contents
    if (!hasError)
    {
        forAll(zones, zi)
        {
            if (zones[zi].checkParallelSync(false))
            {
                hasError = true;

                if (debug || (report && Pstream::master()))
                {
                    Info<< " ***Zone " << zones[zi].name()
                        << " of type " << zones[zi].type()
                        << " is not correctly synchronised"
                        << " across coupled boundaries."
                        << " (coupled faces are either not both"
                        << " present in set or have same flipmap)" << endl;
                }
            }
        }
    }

    return hasError;
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::insert
(
    const List<labelHashSet>& zonesIndices
)
{
    PtrListDictionary<ZoneType>& zones = *this;

    if (zonesIndices.size() != zones.size())
    {
        FatalErrorInFunction
            << "zonesIndices.size() " << zonesIndices.size()
            << " != number of zones " << zones.size()
            << exit(FatalError);
    }

    forAll(zonesIndices, zonei)
    {
        zones[zonei].insert(zonesIndices[zonei]);
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::movePoints
(
    const pointField& p
)
{
    PtrListDictionary<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].movePoints(p);
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::topoChange
(
    const polyTopoChangeMap& map
)
{
    PtrListDictionary<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].topoChange(map);
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::mapMesh
(
    const polyMeshMap& map
)
{
    PtrListDictionary<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].mapMesh(map);
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::distribute
(
    const polyDistributionMap& map
)
{
    PtrListDictionary<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].distribute(map);
    }
}


template<class ZoneType, class ZonesType, class MeshType>
void Foam::ZoneList<ZoneType, ZonesType, MeshType>::swap(ZonesType& otherZones)
{
    clearAddressing();
    otherZones.clearAddressing();

    PtrListDictionary<ZoneType>& zones = *this;

    DynamicList<label> toOtherZone;

    forAll(zones, zi)
    {
        const label ozi = otherZones.findIndex(zones[zi].name());

        if (ozi < 0)
        {
            toOtherZone.append(zi);
        }
    }

    forAll(otherZones, ozi)
    {
        const label zi = this->findIndex(otherZones[ozi].name());

        if (zi < 0)
        {
            zones.append
            (
                otherZones[ozi].name(),
                otherZones[ozi].clone
                (
                    static_cast<const ZonesType&>(*this)
                )
            );
            otherZones.set(ozi, otherZones[ozi].name(), nullptr);
        }
        else
        {
            zones[zi].swap(otherZones[ozi]);
        }
    }

    forAll(toOtherZone, i)
    {
        otherZones.PtrListDictionary<ZoneType>::append
        (
            zones[toOtherZone[i]].name(),
            zones[toOtherZone[i]].clone(otherZones)
        );
        zones.set(toOtherZone[i], zones[toOtherZone[i]].name(), nullptr);
    }

    zones.shrink();
    otherZones.shrink();
}


template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::readIfPresent()
{
    readOpt() = IOobject::READ_IF_PRESENT;
    return read();
}


template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


template<class ZoneType, class ZonesType, class MeshType>
bool Foam::ZoneList<ZoneType, ZonesType, MeshType>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    if (this->size())
    {
        return regIOobject::writeObject(fmt, ver, cmp, write);
    }
    else
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ZoneType, class ZonesType, class MeshType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ZoneList<ZoneType, ZonesType, MeshType>& zones
)
{
    os  << zones.size() << nl << token::BEGIN_LIST;

    forAll(zones, zi)
    {
        zones[zi].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// ************************************************************************* //
