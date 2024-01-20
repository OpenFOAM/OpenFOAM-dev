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
    Foam::functionObjectList

Description
    List of function objects with start(), execute() and end() functions
    that is called for each object.

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    functionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectList_H
#define functionObjectList_H

#include "IOdictionary.H"
#include "PtrList.H"
#include "functionObject.H"
#include "SHA1Digest.H"
#include "includeFuncEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;
class argList;

/*---------------------------------------------------------------------------*\
                     Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    public IOdictionary,
    private PtrList<functionObject>
{
    // Private Data

        //- A list of SHA1 digests for the function object dictionaries
        List<SHA1Digest> digests_;

        //- Quick lookup of the index into functions/digests
        HashTable<label> indices_;

        const Time& time_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        //- Tracks if read() was called while execution is on
        bool updated_;


    // Private Member Functions

        //- Read the function dictionary if present
        //  otherwise read the optional functions entry in controlDict
        typeIOobject<IOdictionary> readFunctionsDict
        (
            const Time& t,
            const bool execution
        );

        //- Read the functions dictionary and set construct the functionObjects
        virtual bool readDict();

        //- Remove and return the function object pointer by name,
        //  and returns the old index via the parameter.
        //  Returns a nullptr (and index -1) if it didn't exist
        functionObject* remove(const word&, label& oldIndex);


public:

    // Constructors

        //- Construct from Time and the execution setting.
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time& runTime,
            const bool execution=true
        );

        //- Disallow default bitwise copy construction
        functionObjectList(const functionObjectList&) = delete;

        //- Construct and return a functionObjectList for an application.
        //  If the "dict" argument is specified the functionObjectList is
        //  constructed from that dictionary which is returned as
        //  controlDict otherwise the functionObjectList is constructed
        //  from the "functions" sub-dictionary of "system/controlDict"
        static autoPtr<functionObjectList> New
        (
            const argList& args,
            const Time& runTime
        );


    //- Destructor
    virtual ~functionObjectList();


    // Member Functions

        //- Return the number of elements in the List.
        using PtrList<functionObject>::size;

        //- Return true if the List is empty (ie, size() is zero).
        using PtrList<functionObject>::empty;

        //- Access to the functionObjects
        using PtrList<functionObject>::operator[];

        //- Allow addition of functionObjects to list
        using PtrList<functionObject>::append;

        //- Clear the list of function objects
        void clear();

        //- Find the ID of a given function object by name
        label findObjectID(const word& name) const;

        //- Read and set the function objects if their data have changed
        virtual bool read();

        //- Called at the start of the time-loop
        bool start();

        //- Called at each ++ or += of the time-loop.
        //  foamPostProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        bool execute();

        //- Called when Time::run() determines that the time-loop exits
        bool end();

        //- Return the time to the next write
        scalar timeToNextAction();

        //- Return the maximum time-step for stable operation
        virtual scalar maxDeltaT() const;

        //- Update topology using the given map
        virtual void movePoints(const polyMesh& mesh);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const functionObjectList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
