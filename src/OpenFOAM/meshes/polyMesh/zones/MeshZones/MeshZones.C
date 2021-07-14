/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "entry.H"
#include "demandDrivenData.H"
#include "stringListOps.H"
#include "Pstream.H"

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

        forAll(*this, zoneI)
        {
            nObjects += this->operator[](zoneI).size();
        }

        zoneMapPtr_ = new Map<label>(2*nObjects);
        Map<label>& zm = *zoneMapPtr_;

        // Fill in objects of all zones into the map.  The key is the global
        // object index and the result is the zone index
        forAll(*this, zoneI)
        {
            const labelList& zoneObjects = this->operator[](zoneI);

            forAll(zoneObjects, objI)
            {
                zm.insert(zoneObjects[objI], zoneI);
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

        forAll(zones, zoneI)
        {
            zones.set
            (
                zoneI,
                ZoneType::New
                (
                    patchEntries[zoneI].keyword(),
                    patchEntries[zoneI].dict(),
                    zoneI,
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
        forAll(zones, zoneI)
        {
            zones.set(zoneI, mpz[zoneI].clone(*this).ptr());
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

    forAll(zones, zoneI)
    {
        lst[zoneI] = zones[zoneI].type();
    }

    return lst;
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::MeshZones<ZoneType, MeshType>::names() const
{
    const PtrList<ZoneType>& zones = *this;

    wordList lst(zones.size());

    forAll(zones, zoneI)
    {
        lst[zoneI] = zones[zoneI].name();
    }

    return lst;
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
Foam::label Foam::MeshZones<ZoneType, MeshType>::findIndex
(
    const wordRe& key
) const
{
    if (!key.empty())
    {
        if (key.isPattern())
        {
            labelList indices = this->findIndices(key);

            // return first element
            if (!indices.empty())
            {
                return indices[0];
            }
        }
        else
        {
            forAll(*this, i)
            {
                if (key == operator[](i).name())
                {
                    return i;
                }
            }
        }
    }

    // not found
    return -1;
}


template<class ZoneType, class MeshType>
Foam::label Foam::MeshZones<ZoneType, MeshType>::findZoneID
(
    const word& zoneName
) const
{
    const PtrList<ZoneType>& zones = *this;

    forAll(zones, zoneI)
    {
        if (zones[zoneI].name() == zoneName)
        {
            return zoneI;
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
void Foam::MeshZones<ZoneType, MeshType>::clearAddressing()
{
    deleteDemandDrivenData(zoneMapPtr_);

    PtrList<ZoneType>& zones = *this;

    forAll(zones, zoneI)
    {
        zones[zoneI].clearAddressing();
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

    forAll(zones, zoneI)
    {
        inError |= zones[zoneI].checkDefinition(report);
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
        forAll(zones, zoneI)
        {
            if (zones[zoneI].checkParallelSync(false))
            {
                hasError = true;

                if (debug || (report && Pstream::master()))
                {
                    Info<< " ***Zone " << zones[zoneI].name()
                        << " of type " << zones[zoneI].type()
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

    forAll(zones, zoneI)
    {
        zones[zoneI].movePoints(p);
    }
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
    const label zoneI = findZoneID(zoneName);

    if (zoneI < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zoneI);
}


template<class ZoneType, class MeshType>
ZoneType& Foam::MeshZones<ZoneType, MeshType>::operator[]
(
    const word& zoneName
)
{
    const label zoneI = findZoneID(zoneName);

    if (zoneI < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zoneI);
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

    forAll(zones, zoneI)
    {
        zones[zoneI].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// ************************************************************************* //
