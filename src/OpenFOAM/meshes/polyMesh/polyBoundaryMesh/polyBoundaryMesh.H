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

        mutable autoPtr<labelList> patchIndicesPtr_;

        mutable autoPtr<labelList> patchFaceIndicesPtr_;

        mutable autoPtr<HashTable<labelList, word>> groupPatchIndicesPtr_;

        //- Edges of neighbouring patches
        mutable autoPtr<List<labelPairList>> nbrEdgesPtr_;


    // Private Member Functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();


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

        //- Move constructor
        polyBoundaryMesh(polyBoundaryMesh&&) = default;

        //- Disallow default bitwise copy construction
        polyBoundaryMesh(const polyBoundaryMesh&) = delete;


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
        const List<labelPairList>& nbrEdges() const;

        //- Return the list of patch names
        wordList toc() const;

        //- Return the list of patch names
        wordList names() const;

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of physical types
        wordList physicalTypes() const;

        //- Return patch indices for all matches. Optionally matches patchGroups
        labelList findIndices
        (
            const wordRe&,
            const bool usePatchGroups = true
        ) const;

        //- Find patch index given a name
        label findIndex(const word& patchName) const;

        //- Find patch indices for a given polyPatch type
        template<class Type>
        labelHashSet findIndices() const;

        //- Return patch index for a given face label
        label whichPatch(const label faceIndex) const;

        //- Boundary face patch indices
        const labelList& patchIndices() const;

        //- Boundary face patch face indices
        const labelList& patchFaceIndices() const;

        //- Per patch group the patch indices
        const HashTable<labelList, word>& groupPatchIndices() const;

        //- Set/add group with patches
        void setGroup(const word& groupName, const labelList& patchIndices);

        //- Return the set of patch indices corresponding to the given names
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
            const labelUList& patchIndices,
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
        void topoChange();

        //- Rename the patches. If validBoundary is set this calls topoChange()
        //  after reordering to recalculate data (so the call needs to be
        //  parallel synchronised).
        void renamePatches(const wordUList& newNames, const bool validBoundary);

        //- Reorders the patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique (is a
        //  shuffle). If validBoundary is set this calls topoChange() after
        //  reordering to recalculate data (so the call needs to be parallel
        //  synchronised).
        void reorderPatches
        (
            const labelUList& newToOld,
            const bool validBoundary
        );

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Write using given format, version and form uncompression
        bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;


    // Member Operators

        //- Return const and non-const reference to polyPatch by index.
        using polyPatchList::operator[];

        //- Return const reference to polyPatch by name.
        const polyPatch& operator[](const word&) const;

        //- Return reference to polyPatch by name.
        polyPatch& operator[](const word&);

        //- Disallow default bitwise assignment
        void operator=(const polyBoundaryMesh&) = delete;


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
