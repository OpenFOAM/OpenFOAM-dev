/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Global
    postProcess

Description
    Execute application functionObjects to post-process existing results.

    If the "dict" argument is specified the functionObjectList is constructed
    from that dictionary otherwise the functionObjectList is constructed from
    the "functions" sub-dictionary of "system/controlDict"

    Multiple time-steps may be processed and the standard utility time
    controls are provided.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef CREATE_TIME
    #define CREATE_TIME createTime.H
#endif

#ifndef CREATE_MESH
    #define CREATE_MESH createMesh.H
#endif

#ifndef CREATE_FIELDS
    #define CREATE_FIELDS createFields.H
#endif

#ifndef CREATE_CONTROL
    #define CREATE_CONTROL createControl.H
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define INCLUDE_FILE(X) INCLUDE_FILE2(X)
#define INCLUDE_FILE2(X) #X

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::argList::addBoolOption
(
    argList::postProcessOptionName,
    "Execute functionObjects only"
);

if (argList::postProcess(argc, argv))
{
    Foam::timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "addDictOption.H"

    #include "setRootCase.H"
    #include INCLUDE_FILE(CREATE_TIME)
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include INCLUDE_FILE(CREATE_MESH)

    #ifndef NO_CONTROL
    #include INCLUDE_FILE(CREATE_CONTROL)
    #endif

    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary functionObjectsDict;

    // Construct functionObjectList
    autoPtr<functionObjectList> functionObjectsPtr
    (
        functionObjectList::New(args, runTime, functionObjectsDict)
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        if (mesh.readUpdate() != polyMesh::UNCHANGED)
        {
            // Update functionObjects if mesh changes
            functionObjectsPtr =
                functionObjectList::New(args, runTime, functionObjectsDict);
        }

        #include INCLUDE_FILE(CREATE_FIELDS)

        #ifdef CREATE_FIELDS_2
        #include INCLUDE_FILE(CREATE_FIELDS_2)
        #endif

        #ifdef CREATE_FIELDS_3
        #include INCLUDE_FILE(CREATE_FIELDS_3)
        #endif

        FatalIOError.throwExceptions();

        try
        {
            functionObjectsPtr->execute(true);
        }
        catch (IOerror& err)
        {
            Warning<< err << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef INCLUDE_FILE
#undef INCLUDE_FILE2

#undef CREATE_MESH
#undef CREATE_FIELDS
#undef CREATE_CONTROL

// ************************************************************************* //
