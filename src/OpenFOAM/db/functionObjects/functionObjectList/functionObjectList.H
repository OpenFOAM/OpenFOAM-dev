/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "PtrList.H"
#include "functionObject.H"
#include "SHA1Digest.H"
#include "HashTable.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class argList;

/*---------------------------------------------------------------------------*\
                     Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    private PtrList<functionObject>
{
    // Private data

        //- A list of SHA1 digests for the function object dictionaries
        List<SHA1Digest> digests_;

        //- Quick lookup of the index into functions/digests
        HashTable<label> indices_;

        const Time& time_;

        //- The parent dictionary containing a "functions" entry
        //  This entry can either be a list or a dictionary of
        //  functionObject specifications
        const dictionary& parentDict_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        //- Tracks if read() was called while execution is on
        bool updated_;


    // Private Member Functions

        //- Remove and return the function object pointer by name,
        //  and returns the old index via the parameter.
        //  Returns a nullptr (and index -1) if it didn't exist
        functionObject* remove(const word&, label& oldIndex);

        //- Search the specified directory for functionObject
        //  configuration files, add to the given map and recurse
        static void listDir(const fileName& dir, HashSet<word>& foMap);

        //- Disallow default bitwise copy construct
        functionObjectList(const functionObjectList&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectList&);


public:

    // Static data members

        //- Default relative path to the directory structure
        //  containing the functionObject dictionary files
        static fileName functionObjectDictPath;


    // Constructors

        //- Construct from Time and the execution setting.
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time& runTime,
            const bool execution=true
        );

        //- Construct from Time, a dictionary with "functions" entry
        //  and the execution setting.
        //  \param[in]  t - the other Time instance to construct from
        //  \param[in]  parentDict - the parent dictionary containing
        //    a "functions" entry, which can either be a list or a dictionary
        //    of functionObject specifications.
        //  \param[in]  execution - whether the function objects should execute
        //    or not. Default: true.
        functionObjectList
        (
            const Time& runTime,
            const dictionary& parentDict,
            const bool execution=true
        );

        //- Construct and return a functionObjectList for an application.
        //  If the "dict" argument is specified the functionObjectList is
        //  constructed from that dictionary which is returned as
        //  controlDict otherwise the functionObjectList is constructed
        //  from the "functions" sub-dictionary of "system/controlDict"
        static autoPtr<functionObjectList> New
        (
            const argList& args,
            const Time& runTime,
            dictionary& controlDict,
            HashSet<word>& requiredFields
        );


    //- Destructor
    ~functionObjectList();


    // Member Functions

        //- Return the number of elements in the List.
        using PtrList<functionObject>::size;

        //- Return true if the List is empty (ie, size() is zero).
        using PtrList<functionObject>::empty;

        //- Access to the functionObjects
        using PtrList<functionObject>::operator[];

        //- Clear the list of function objects
        void clear();

        //- Find the ID of a given function object by name
        label findObjectID(const word& name) const;

        //- Print a list of functionObject configuration files in
        //  user/group/shipped directories.
        //  The search scheme allows for version-specific and
        //  version-independent files using the following hierarchy:
        //  - \b user settings:
        //    - ~/.OpenFOAM/\<VERSION\>/caseDicts/postProcessing
        //    - ~/.OpenFOAM/caseDicts/postProcessing
        //  - \b group (site) settings (when $WM_PROJECT_SITE is set):
        //    - $WM_PROJECT_SITE/\<VERSION\>/caseDicts/postProcessing
        //    - $WM_PROJECT_SITE/caseDicts/postProcessing
        //  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
        //    - $WM_PROJECT_INST_DIR/site/\<VERSION\>/caseDicts/postProcessing
        //    - $WM_PROJECT_INST_DIR/site/caseDicts/postProcessing
        //  - \b other (shipped) settings:
        //    - $WM_PROJECT_DIR/etc/caseDicts/postProcessing
        static void list();

        //- Search for functionObject dictionary file for given region
        //  and the user/group/shipped directories.
        //  The search scheme allows for version-specific and
        //  version-independent files using the following hierarchy:
        //  - \b user settings:
        //    - ~/.OpenFOAM/\<VERSION\>/caseDicts/postProcessing
        //    - ~/.OpenFOAM/caseDicts/postProcessing
        //  - \b group (site) settings (when $WM_PROJECT_SITE is set):
        //    - $WM_PROJECT_SITE/\<VERSION\>/caseDicts/postProcessing
        //    - $WM_PROJECT_SITE/caseDicts/postProcessing
        //  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
        //    - $WM_PROJECT_INST_DIR/site/\<VERSION\>/caseDicts/postProcessing
        //    - $WM_PROJECT_INST_DIR/site/caseDicts/postProcessing
        //  - \b other (shipped) settings:
        //    - $WM_PROJECT_DIR/etc/caseDicts/postProcessing
        //
        //  \return The path of the functionObject dictionary file if found
        //  otherwise null
        static fileName findRegionDict
        (
            const word& funcPath,
            const word& region = word::null
        );

        //- Search for functionObject dictionary file for given region
        //  and if not present also search the case directory as well as the
        //  user/group/shipped directories.
        //  The search scheme allows for version-specific and
        //  version-independent files using the following hierarchy:
        //  - \b user settings:
        //    - ~/.OpenFOAM/\<VERSION\>/caseDicts/postProcessing
        //    - ~/.OpenFOAM/caseDicts/postProcessing
        //  - \b group (site) settings (when $WM_PROJECT_SITE is set):
        //    - $WM_PROJECT_SITE/\<VERSION\>/caseDicts/postProcessing
        //    - $WM_PROJECT_SITE/caseDicts/postProcessing
        //  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
        //    - $WM_PROJECT_INST_DIR/site/\<VERSION\>/caseDicts/postProcessing
        //    - $WM_PROJECT_INST_DIR/site/caseDicts/postProcessing
        //  - \b other (shipped) settings:
        //    - $WM_PROJECT_DIR/etc/caseDicts/postProcessing
        //
        //  \return The path of the functionObject dictionary file if found
        //  otherwise null
        static fileName findDict
        (
            const word& funcName,
            const word& region = word::null
        );

        //- Read the specified functionObject configuration dictionary parsing
        //  the optional arguments included in the name 'funcNameArgs0',
        //  inserting 'field' or 'fields' entries as required and merging the
        //  resulting functionObject dictionary into 'functionsDict'.  Any
        //  fields required to execute the functionObject are added to
        //  'requiredFields'
        static bool readFunctionObject
        (
            const string& funcNameArgs0,
            dictionary& functionsDict,
            HashSet<word>& requiredFields,
            const word& region = word::null
        );

        //- Read and set the function objects if their data have changed
        bool read();

        //- Switch the function objects on
        void on();

        //- Switch the function objects off
        void off();

        //- Return the execution status (on/off) of the function objects
        bool status() const;

        //- Called at the start of the time-loop
        bool start();

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        bool execute();

        //- Called when Time::run() determines that the time-loop exits
        bool end();

        //- Override the time-step value
        bool setTimeStep();

        //- Return the time to the next write
        scalar timeToNextWrite();

        //- Update for changes of mesh
        void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
