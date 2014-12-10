/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::polyBoundaryMesh

Description
    Foam::polyBoundaryMesh

SourceFiles
    polyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMesh_H
#define polyBoundaryMesh_H

#include "polyPatchList.H"
#include "regIOobject.H"
#include "labelPair.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class wordRe;

// Forward declaration of friend functions and operators

Ostream& operator<<(Ostream&, const polyBoundaryMesh&);


/*---------------------------------------------------------------------------*\
                      Class polyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class polyBoundaryMesh
:
    public polyPatchList,
    public regIOobject
{
    // private data

        //- Reference to mesh
        const polyMesh& mesh_;

        mutable autoPtr<labelList> patchIDPtr_;

        mutable autoPtr<HashTable<labelList, word> > groupPatchIDsPtr_;

        //- Edges of neighbouring patches
        mutable autoPtr<List<labelPairList> > neighbourEdgesPtr_;


    // Private Member Functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- Disallow construct as copy
        polyBoundaryMesh(const polyBoundaryMesh&);

        //- Disallow assignment
        void operator=(const polyBoundaryMesh&);


public:

    //- Declare friendship with polyMesh
    friend class polyMesh;


    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        //- Read constructor given IOobject and a polyMesh reference
        //  Note point pointers are unset, only used in copying meshes
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&
        );

        //- Construct given size
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&,
            const label size
        );

        //- Construct given polyPatchList
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&,
            const polyPatchList&
        );


    //- Destructor
    ~polyBoundaryMesh();

        //- Clear geometry at this level and at patches
        void clearGeom();

        //- Clear addressing at this level and at patches
        void clearAddressing();


    // Member Functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Per patch the edges on the neighbouring patch. Is for every external
        //  edge the neighbouring patch and neighbouring (external) patch edge
        //  label. Note that edge indices are offset by nInternalEdges to keep
        //  it as much as possible consistent with coupled patch addressing
        //  (where coupling is by local patch face index).
        //  Only valid for singly connected polyBoundaryMesh and not parallel
        const List<labelPairList>& neighbourEdges() const;

        //- Return a list of patch names
        wordList names() const;

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of physical types
        wordList physicalTypes() const;

        //- Return patch indices for all matches. Optionally matches patchGroups
        labelList findIndices
        (
            const keyType&,
            const bool usePatchGroups = true
        ) const;

        //- Return patch index for the first match, return -1 if not found
        label findIndex(const keyType&) const;

        //- Find patch index given a name
        label findPatchID(const word& patchName) const;

        //- Find patch indices for a given polyPatch type
        template<class Type>
        labelHashSet findPatchIDs() const;

        //- Return patch index for a given face label
        label whichPatch(const label faceIndex) const;

        //- Per boundary face label the patch index
        const labelList& patchID() const;

        //- Per patch group the patch indices
        const HashTable<labelList, word>& groupPatchIDs() const;

        //- Set/add group with patches
        void setGroup(const word& groupName, const labelList& patchIDs);

        //- Return the set of patch IDs corresponding to the given names
        //  By default warns if given names are not found. Optionally
        //  matches to patchGroups as well as patchNames
        labelHashSet patchSet
        (
            const UList<wordRe>& patchNames,
            const bool warnNotFound = true,
            const bool usePatchGroups = true
        ) const;

        //- Match the patches to groups. Returns all the (fully matched) groups
        //  and any remaining unmatched patches.
        void matchGroups
        (
            const labelUList& patchIDs,
            wordList& groups,
            labelHashSet& nonGroupPatches
        ) const;

        //- Check whether all procs have all patches and in same order. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Check boundary definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Correct polyBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- Correct polyBoundaryMesh after topology update
        void updateMesh();

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle) If validBoundary calls updateMesh()
        //  after reordering to recalculate data (so call needs to be parallel
        //  sync in that case)
        void reorder(const labelUList&, const bool validBoundary);

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Write using given format, version and form uncompression
        bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;


    // Member Operators

        //- Return const and non-const reference to polyPatch by index.
        using polyPatchList::operator[];

        //- Return const reference to polyPatch by name.
        const polyPatch& operator[](const word&) const;

        //- Return reference to polyPatch by name.
        polyPatch& operator[](const word&);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "polyBoundaryMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
