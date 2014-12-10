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

Class
    Foam::topoSet

Description
    General set of labels of mesh quantity (points, cells, faces).

    Contains various 'notImplemented' functions, but I do not want to make
    this class abstract since it is quite handy to work on topoSets.

SourceFiles
    topoSet.C

\*---------------------------------------------------------------------------*/

#ifndef topoSet_H
#define topoSet_H

#include "HashSet.H"
#include "regIOobject.H"
#include "labelList.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class mapPolyMesh;
class polyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class topoSet Declaration
\*---------------------------------------------------------------------------*/

class topoSet
:
    public regIOobject,
    public labelHashSet
{

protected:

    // Protected Member Functions

        //- Update map from map. Used to update cell/face labels
        //  after morphing
        void updateLabels(const labelList& map);

        //- Check validity of contents.
        void check(const label maxLabel);

        //- Write part of contents nicely formatted. Prints labels only.
        void writeDebug
        (
            Ostream& os,
            const label maxElem,
            topoSet::const_iterator& iter,
            label& elemI
        ) const;

        //- Write part of contents nicely formatted. Prints label
        //  and corresponding coordinate.
        void writeDebug
        (
            Ostream& os,
            const pointField& coords,
            const label maxElem,
            topoSet::const_iterator& iter,
            label& elemI
        ) const;

        //- Write labels and coordinates columnwise to os. Truncate to maxLen.
        void writeDebug
        (
            Ostream& os,
            const pointField& coords,
            const label maxLen
        ) const;


        //- Disallow default bitwise copy construct
        topoSet(const topoSet&);

public:

    //- Runtime type information
    TypeName("topoSet");


    // Static

        //- Name of file set will use.
        static fileName localPath(const polyMesh& mesh, const word& name);


    // Declare run-time constructor selection table

        // For the direct constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSet,
            word,
            (
                const polyMesh& mesh,
                const word& name,
                readOption r,
                writeOption w
            ),
            (mesh, name, r, w)
        );

        // For the constructor from size
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSet,
            size,
            (
                const polyMesh& mesh,
                const word& name,
                const label size,
                writeOption w
            ),
            (mesh, name, size, w)
        );

        // For the constructor as copy
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSet,
            set,
            (
                const polyMesh& mesh,
                const word& name,
                const topoSet& set,
                writeOption w
            ),
            (mesh, name, set, w)
        );


    // Constructors


        //- Construct from IOobject as explicitly passed type.
        //  Can't use typeName info here since subclasses not yet instantiated
        topoSet(const IOobject&, const word& wantedType);

        //- Construct from polyMesh and name. Searches for a polyMesh/sets
        //  directory but not beyond mesh.facesInstance.
        topoSet
        (
            const polyMesh& mesh,
            const word& wantedType,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct empty from additional size of labelHashSet.
        //  Searches for a polyMesh/sets
        //  directory but not beyond mesh.facesInstance.
        topoSet
        (
            const polyMesh& mesh,
            const word& name,
            const label,
            writeOption w=NO_WRITE
        );

        //- Construct empty from additional labelHashSet
        //  Searches for a polyMesh/sets
        //  directory but not beyond mesh.facesInstance.
        topoSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet&,
            writeOption w=NO_WRITE
        );

        //- Construct empty from IOobject and size.
        topoSet(const IOobject&, const label size);

        //- Construct from IOobject and labelHashSet.
        topoSet(const IOobject&, const labelHashSet&);



        //- Clone
        autoPtr<topoSet> clone() const
        {
            notImplemented("autoPtr<topoSet> clone() const");
            return autoPtr<topoSet>(NULL);
        }


    // Selectors

        //- Return a pointer to a toposet read from file
        static autoPtr<topoSet> New
        (
            const word& setType,
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Return a pointer to a new toposet of given size
        static autoPtr<topoSet> New
        (
            const word& setType,
            const polyMesh& mesh,
            const word& name,
            const label size,
            writeOption w=NO_WRITE
        );

        //- Return a pointer to a new toposet as copy of another toposet
        static autoPtr<topoSet> New
        (
            const word& setType,
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            writeOption w=NO_WRITE
        );


    //- Destructor
    virtual ~topoSet();


    // Member functions

        //- Invert contents. (insert all members 0..maxLen-1 which were not in
        //  set)
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add elements present in set.
        virtual void addSet(const topoSet& set);

        //- Delete elements present in set.
        virtual void deleteSet(const topoSet& set);

        //- Sync set across coupled patches.
        virtual void sync(const polyMesh& mesh);

        //- Write labels columnwise to os. Truncate to maxLen.
        virtual void writeDebug(Ostream& os, const label maxLen) const;

        //- Like above but also writes mesh related quantity
        //  (usually coordinate).
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const = 0;

        //- Write contents.
        virtual bool writeData(Ostream&) const;

        //- Update any stored data for new labels. Not implemented.
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Return max allowable index (+1). Not implemented.
        virtual label maxSize(const polyMesh& mesh) const = 0;



    // Member operators

        //- Copy labelHashSet part only
        void operator=(const topoSet&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
