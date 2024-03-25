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
    Foam::MeshZones

Description
    A list of mesh zones.

SourceFiles
    MeshZones.C

\*---------------------------------------------------------------------------*/

#ifndef MeshZones_H
#define MeshZones_H

#include "regIOobject.H"
#include "pointField.H"
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

template<class ZoneType, class MeshType> class MeshZones;

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream&, const MeshZones<ZoneType, MeshType>&);

/*---------------------------------------------------------------------------*\
                          Class MeshZones Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType, class MeshType>
class MeshZones
:
    public PtrList<ZoneType>,
    public regIOobject
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
        MeshZones
        (
            const IOobject&,
            const MeshType&
        );

        //- Construct given size
        MeshZones
        (
            const IOobject&,
            const MeshType&,
            const label size
        );

         //- Construct given a PtrList
        MeshZones
        (
            const IOobject&,
            const MeshType&,
            const PtrList<ZoneType>&
        );

        //- Move constructor
        MeshZones(MeshZones&&) = default;

        //- Disallow default bitwise copy construction
        MeshZones(const MeshZones&) = delete;


    //- Destructor
    ~MeshZones();


    // Member Functions

        //- Return the mesh reference
        const MeshType& mesh() const
        {
            return mesh_;
        }

        //- Return true if objectIndex is in any zone
        bool found(const label objectIndex) const;

        //- Given a global object index, return the list of zones it is in
        labelList whichZones(const label objectIndex) const;

        //- Return the list of flips for the facei
        //  for each of the face zones in faceiZones
        boolList zonesFlipFace
        (
            const label facei,
            const labelList& faceiZones
        );

        //- Return a list of zone types
        wordList types() const;

        //- Return a list of zone names
        wordList names() const;

        //- Return true if the given zoneName is present
        bool found(const word& zoneName) const;

        //- Find the zone index given the zone name
        label findIndex(const word& zoneName) const;

        //- Find and return the zone indices for all matches
        labelList findIndices(const wordRe&) const;

        //- Mark cells that match the zone specification
        PackedBoolList findMatching(const wordRe&) const;

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
        void swap(MeshZones&);

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

        //- Return const and non-const reference to ZoneType by index.
        using PtrList<ZoneType>::operator[];

        //- Return const reference to ZoneType by name.
        const ZoneType& operator[](const word&) const;

        //- Return reference to ZoneType by name.
        ZoneType& operator[](const word&);

        //- Disallow default bitwise assignment
        void operator=(const MeshZones<ZoneType, MeshType>&) = delete;


    // Ostream operator

        friend Ostream& operator<< <ZoneType, MeshType>
        (
            Ostream&,
            const MeshZones<ZoneType, MeshType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshZones.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
