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

Class
    Foam::ZoneList

Description
    A list of mesh zones.

SourceFiles
    ZoneList.C

\*---------------------------------------------------------------------------*/

#ifndef ZoneList_H
#define ZoneList_H

#include "PtrListDictionary.H"
#include "regIOobject.H"
#include "pointFieldFwd.H"
#include "Map.H"
#include "boolList.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

template<class ZoneType, class ZonesType, class MeshType> class ZoneList;

template<class ZoneType, class ZonesType, class MeshType>
Ostream& operator<<(Ostream&, const ZoneList<ZoneType, ZonesType, MeshType>&);

/*---------------------------------------------------------------------------*\
                          Class ZoneList Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType, class ZonesType, class MeshType>
class ZoneList
:
    public regIOobject,
    public PtrListDictionary<ZoneType>
{
    // Private Data

        //- Reference to mesh
        const MeshType& mesh_;


    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool read();


public:

    // Constructors

        //- Read constructor given IOobject and a MeshType reference
        ZoneList
        (
            const IOobject&,
            const MeshType&
        );

        //- Move constructor
        ZoneList(ZoneList&&) = default;

        //- Disallow default bitwise copy construction
        ZoneList(const ZoneList&) = delete;


    //- Destructor
    ~ZoneList();


    // Member Functions

        //- Return the mesh reference
        const MeshType& mesh() const
        {
            return mesh_;
        }

        using PtrListDictionary<ZoneType>::found;

        //- Return true if objectIndex is in any zone
        bool found(const label objectIndex) const;

        //- Given a global object index, return the list of zones it is in
        labelList whichZones(const label objectIndex) const;

        //- Return a list of zone types
        wordList types() const;

        //- Return the set of zone indices corresponding to the given names
        //  By default warns if given names are not found
        labelHashSet zoneSet
        (
            const UList<wordRe>& zoneNames,
            const bool warnNotFound = true
        ) const;

        //- Append or update a zone
        void append(ZoneType*) const;

        //- Append or update a zone
        void append(const ZoneType&) const;

        //- Clear addressing
        void clearAddressing();

        //- Clear the zones
        void clear();

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have all zones and in same order. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Insert given indices into zones
        void insert(const List<labelHashSet>& zonesIndices);

        //- Correct zones after moving points
        virtual void movePoints(const pointField&);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Swap zones
        //  For run-time mesh replacement and mesh to mesh mapping
        void swap(ZonesType&);

        //- Read zones if the zones file is present
        bool readIfPresent();

        //- writeData member function required by regIOobject
        virtual bool writeData(Ostream&) const;

        //- Write using given format, version and compression
        //  Only write the zones file if there are zones in the list
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write
        ) const;


    // Member Operators

        //- Return const and non-const reference to ZoneType by index.
        using PtrListDictionary<ZoneType>::operator[];

        //- Disallow default bitwise assignment
        void operator=(const ZoneList<ZoneType, ZonesType, MeshType>&) = delete;


    // Ostream operator

        friend Ostream& operator<< <ZoneType, ZonesType, MeshType>
        (
            Ostream&,
            const ZoneList<ZoneType, ZonesType, MeshType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ZoneList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
