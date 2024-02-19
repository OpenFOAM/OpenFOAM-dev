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

#include "MeshZones.H"
#include "Pstream.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::calcZoneMap() const
{
    // It is an error to attempt to recalculate cellEdges
    // if the pointer is already set
    if (zoneMapPtr_)
    {
        FatalErrorInFunction
            << "zone map already calculated"
            << abort(FatalError);
    }
    else
    {
        // Count number of objects in all zones
        label nObjects = 0;

        forAll(*this, zi)
        {
            nObjects += this->operator[](zi).size();
        }

        zoneMapPtr_ = new Map<label>(2*nObjects);
        Map<label>& zm = *zoneMapPtr_;

        // Fill in objects of all zones into the map.  The key is the global
        // object index and the result is the zone index
        forAll(*this, zi)
        {
            const labelList& zoneObjects = this->operator[](zi);

            forAll(zoneObjects, objI)
            {
                zm.insert(zoneObjects[objI], zi);
            }
        }
    }
}


template<class ZoneType, class MeshType>
bool Foam::MeshZones<ZoneType, MeshType>::read()
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

        PtrList<ZoneType>& zones = *this;

        // Read zones
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        zones.setSize(patchEntries.size());

        forAll(zones, zi)
        {
            zones.set
            (
                zi,
                ZoneType::New
                (
                    patchEntries[zi].keyword(),
                    patchEntries[zi].dict(),
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "MeshZones::MeshZones"
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

template<class ZoneType, class MeshType>
Foam::MeshZones<ZoneType, MeshType>::MeshZones
(
    const IOobject& io,
    const MeshType& mesh
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(nullptr)
{
    read();
}


template<class ZoneType, class MeshType>
Foam::MeshZones<ZoneType, MeshType>::MeshZones
(
    const IOobject& io,
    const MeshType& mesh,
    const label size
)
:
    PtrList<ZoneType>(size),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(nullptr)
{
    // Optionally read contents, otherwise keep size
    read();
}


template<class ZoneType, class MeshType>
Foam::MeshZones<ZoneType, MeshType>::MeshZones
(
    const IOobject& io,
    const MeshType& mesh,
    const PtrList<ZoneType>& mpz
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(nullptr)
{
    if (!read())
    {
        // Nothing read. Use supplied zones
        PtrList<ZoneType>& zones = *this;
        zones.setSize(mpz.size());
        forAll(zones, zi)
        {
            zones.set(zi, mpz[zi].clone(*this).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::MeshZones<ZoneType, MeshType>::~MeshZones()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
const Foam::Map<Foam::label>&
Foam::MeshZones<ZoneType, MeshType>::zoneMap() const
{
    if (!zoneMapPtr_)
    {
        calcZoneMap();
    }

    return *zoneMapPtr_;
}


template<class ZoneType, class MeshType>
Foam::label Foam::MeshZones<ZoneType, MeshType>::whichZone
(
    const label objectIndex
) const
{
    const Map<label>& zm = zoneMap();
    Map<label>::const_iterator zmIter = zm.find(objectIndex);

    if (zmIter == zm.end())
    {
        return -1;
    }
    else
    {
        return zmIter();
    }
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::MeshZones<ZoneType, MeshType>::types() const
{
    const PtrList<ZoneType>& zones = *this;

    wordList lst(zones.size());

    forAll(zones, zi)
    {
        lst[zi] = zones[zi].type();
    }

    return lst;
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::MeshZones<ZoneType, MeshType>::names() const
{
    const PtrList<ZoneType>& zones = *this;

    wordList lst(zones.size());

    forAll(zones, zi)
    {
        lst[zi] = zones[zi].name();
    }

    return lst;
}


template<class ZoneType, class MeshType>
bool Foam::MeshZones<ZoneType, MeshType>::found
(
    const word& zoneName
) const
{
    if (zoneName != word::null)
    {
        forAll(*this, i)
        {
            if (zoneName == operator[](i).name())
            {
                return true;
            }
        }
    }

    // Not found
    return false;
}


template<class ZoneType, class MeshType>
Foam::label Foam::MeshZones<ZoneType, MeshType>::findIndex
(
    const word& zoneName
) const
{
    const PtrList<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        if (zones[zi].name() == zoneName)
        {
            return zi;
        }
    }

    // Zone not found
    if (debug)
    {
        InfoInFunction
            << "Zone named " << zoneName << " not found.  "
            << "List of available zone names: " << names() << endl;
    }

    // not found
    return -1;
}


template<class ZoneType, class MeshType>
Foam::labelList Foam::MeshZones<ZoneType, MeshType>::findIndices
(
    const wordRe& key
) const
{
    labelList indices;

    if (!key.empty())
    {
        if (key.isPattern())
        {
            indices = findStrings(key, this->names());
        }
        else
        {
            indices.setSize(this->size());
            label nFound = 0;
            forAll(*this, i)
            {
                if (key == operator[](i).name())
                {
                    indices[nFound++] = i;
                }
            }
            indices.setSize(nFound);
        }
    }

    return indices;
}


template<class ZoneType, class MeshType>
Foam::PackedBoolList Foam::MeshZones<ZoneType, MeshType>::findMatching
(
    const wordRe& key
) const
{
    PackedBoolList lst;

    const labelList indices = this->findIndices(key);
    forAll(indices, i)
    {
        lst |= static_cast<const labelList&>(this->operator[](indices[i]));
    }

    return lst;
}


template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::append(ZoneType* zonePtr) const
{
    MeshZones<ZoneType, MeshType>& zones =
        const_cast<MeshZones<ZoneType, MeshType>&>(*this);

    if (found(zonePtr->name()))
    {
        zones[zonePtr->name()] = *zonePtr;
        delete zonePtr;
    }
    else
    {
        zones.PtrList<ZoneType>::append(zonePtr);
    }
}


template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::append(const ZoneType& zone) const
{
    MeshZones<ZoneType, MeshType>& zones =
        const_cast<MeshZones<ZoneType, MeshType>&>(*this);

    if (found(zone.name()))
    {
        zones[zone.name()] = zone;
    }
    else
    {
        zones.PtrList<ZoneType>::append(zone.clone(*this));
    }
}


template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::clearAddressing()
{
    deleteDemandDrivenData(zoneMapPtr_);

    PtrList<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].clearAddressing();
    }
}


template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::clear()
{
    clearAddressing();
    PtrList<ZoneType>::clear();
}


template<class ZoneType, class MeshType>
bool Foam::MeshZones<ZoneType, MeshType>::checkDefinition
(
    const bool report
) const
{
    bool inError = false;

    const PtrList<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        inError |= zones[zi].checkDefinition(report);
    }
    return inError;
}


template<class ZoneType, class MeshType>
bool Foam::MeshZones<ZoneType, MeshType>::checkParallelSync
(
    const bool report
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }


    const PtrList<ZoneType>& zones = *this;

    bool hasError = false;

    // Collect all names
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = this->names();
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


template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::movePoints(const pointField& p)
{
    PtrList<ZoneType>& zones = *this;

    forAll(zones, zi)
    {
        zones[zi].movePoints(p);
    }
}


template<class ZoneType, class MeshType>
void Foam::MeshZones<ZoneType, MeshType>::swap(MeshZones& otherZones)
{
    clearAddressing();
    otherZones.clearAddressing();

    PtrList<ZoneType>& zones = *this;

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
        const label zi = findIndex(otherZones[ozi].name());

        if (zi < 0)
        {
            zones.append(otherZones[ozi].clone(*this));
            otherZones.set(ozi, nullptr);
        }
        else
        {
            zones[zi].swap(otherZones[ozi]);
        }
    }

    forAll(toOtherZone, i)
    {
        otherZones.PtrList<ZoneType>::append
        (
            zones[toOtherZone[i]].clone(otherZones)
        );
        zones.set(toOtherZone[i], nullptr);
    }

    zones.shrink();
    otherZones.shrink();
}


template<class ZoneType, class MeshType>
bool Foam::MeshZones<ZoneType, MeshType>::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
const ZoneType& Foam::MeshZones<ZoneType, MeshType>::operator[]
(
    const word& zoneName
) const
{
    const label zi = findIndex(zoneName);

    if (zi < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zi);
}


template<class ZoneType, class MeshType>
ZoneType& Foam::MeshZones<ZoneType, MeshType>::operator[]
(
    const word& zoneName
)
{
    const label zi = findIndex(zoneName);

    if (zi < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zi);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MeshZones<ZoneType, MeshType>& zones
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
