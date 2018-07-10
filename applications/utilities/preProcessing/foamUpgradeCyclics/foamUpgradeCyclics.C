/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Application
    foamUpgradeCyclics

Description
    Tool to upgrade mesh and fields for split cyclics.

Usage
    \b foamUpgradeCyclics [OPTION]

    Options:
      - \par -test
        Suppress writing the updated files with split cyclics

      - \par -enableFunctionEntries
        By default all dictionary preprocessing of fields is disabled

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "entry.H"
#include "IOPtrList.H"
#include "cyclicPolyPatch.H"
#include "dictionaryEntry.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "string.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// Read boundary file without reading mesh
void rewriteBoundary
(
    const bool isTestRun,
    const IOobject& io,
    const fileName& regionPrefix,
    HashTable<word>& thisNames,
    HashTable<word>& nbrNames
)
{
    Info<< "Reading boundary from " << typeFilePath<IOPtrList<entry>>(io)
        << endl;

    // Read PtrList of dictionary.
    const word oldTypeName = IOPtrList<entry>::typeName;
    const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
    IOPtrList<entry> patches(io);
    const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
    // Fake type back to what was in field
    const_cast<word&>(patches.type()) = patches.headerClassName();


    // Replace any 'cyclic'
    label nOldCyclics = 0;
    forAll(patches, patchi)
    {
        const dictionary& patchDict = patches[patchi].dict();

        if (word(patchDict["type"]) == cyclicPolyPatch::typeName)
        {
            if (!patchDict.found("neighbourPatch"))
            {
                Info<< "Patch " << patches[patchi].keyword()
                    << " does not have 'neighbourPatch' entry; assuming it"
                    << " is of the old type." << endl;
                nOldCyclics++;
            }
        }
    }

    Info<< "Detected " << nOldCyclics << " old cyclics." << nl << endl;


    // Save old patches.
    PtrList<entry> oldPatches(patches);

    // Extend
    label nOldPatches = patches.size();
    patches.setSize(nOldPatches+nOldCyclics);

    // Create reordering map
    labelList oldToNew(patches.size());


    // Add new entries
    label addedPatchi = nOldPatches;
    label newPatchi = 0;
    forAll(oldPatches, patchi)
    {
        const dictionary& patchDict = oldPatches[patchi].dict();

        if
        (
            word(patchDict["type"]) == cyclicPolyPatch::typeName
        )
        {
            const word& name = oldPatches[patchi].keyword();

            if (patchDict.found("neighbourPatch"))
            {
                patches.set(patchi, oldPatches.set(patchi, nullptr));
                oldToNew[patchi] = newPatchi++;

                // Check if patches come from automatic conversion
                word oldName;

                string::size_type i = name.rfind("_half0");
                if (i != string::npos)
                {
                    oldName = name.substr(0, i);
                    thisNames.insert(oldName, name);
                    Info<< "Detected converted cyclic patch " << name
                        << " ; assuming it originates from " << oldName
                        << endl;
                }
                else
                {
                    i = name.rfind("_half1");
                    if (i != string::npos)
                    {
                        oldName = name.substr(0, i);
                        nbrNames.insert(oldName, name);
                        Info<< "Detected converted cyclic patch " << name
                            << " ; assuming it originates from " << oldName
                            << endl;
                    }
                }
            }
            else
            {
                label nFaces = readLabel(patchDict["nFaces"]);
                label startFace = readLabel(patchDict["startFace"]);

                Info<< "Detected old style " << word(patchDict["type"])
                    << " patch " << name << " with" << nl
                    << "    nFaces    : " << nFaces << nl
                    << "    startFace : " << startFace << endl;

                word thisName = name + "_half0";
                word nbrName = name + "_half1";

                thisNames.insert(name, thisName);
                nbrNames.insert(name, nbrName);

                // Save current dictionary
                const dictionary patchDict(patches[patchi].dict());

                // Change entry on this side
                patches.set(patchi, oldPatches.set(patchi, nullptr));
                oldToNew[patchi] = newPatchi++;
                dictionary& thisPatchDict = patches[patchi].dict();
                thisPatchDict.add("neighbourPatch", nbrName);
                thisPatchDict.set("nFaces", nFaces/2);
                patches[patchi].keyword() = thisName;

                // Add entry on other side
                patches.set
                (
                    addedPatchi,
                    new dictionaryEntry
                    (
                        nbrName,
                        dictionary::null,
                        patchDict
                    )
                );
                oldToNew[addedPatchi] = newPatchi++;
                dictionary& nbrPatchDict = patches[addedPatchi].dict();
                nbrPatchDict.set("neighbourPatch", thisName);
                nbrPatchDict.set("nFaces", nFaces/2);
                nbrPatchDict.set("startFace", startFace+nFaces/2);
                patches[addedPatchi].keyword() = nbrName;

                Info<< "Replaced with patches" << nl
                    << patches[patchi].keyword() << " with" << nl
                    << "    nFaces    : "
                    << readLabel(thisPatchDict.lookup("nFaces"))
                    << nl
                    << "    startFace : "
                    << readLabel(thisPatchDict.lookup("startFace")) << nl
                    << patches[addedPatchi].keyword() << " with" << nl
                    << "    nFaces    : "
                    << readLabel(nbrPatchDict.lookup("nFaces"))
                    << nl
                    << "    startFace : "
                    << readLabel(nbrPatchDict.lookup("startFace"))
                    << nl << endl;

                addedPatchi++;
            }
        }
        else
        {
            patches.set(patchi, oldPatches.set(patchi, nullptr));
            oldToNew[patchi] = newPatchi++;
        }
    }

    patches.reorder(oldToNew);

    if (returnReduce(nOldCyclics, sumOp<label>()) > 0)
    {
        if (isTestRun)
        {
            // Info<< "-test option: no changes made" << nl << endl;
        }
        else
        {
            if (mvBak(patches.objectPath(), "old"))
            {
                Info<< "Backup to    "
                    << (patches.objectPath() + ".old") << nl;
            }

            Info<< "Write  to    "
                << patches.objectPath() << nl << endl;
            patches.write();
        }
    }
    else
    {
        Info<< "No changes made to boundary file." << nl << endl;
    }
}


void rewriteField
(
    const bool isTestRun,
    const Time& runTime,
    const word& fieldName,
    const HashTable<word>& thisNames,
    const HashTable<word>& nbrNames
)
{
    // Read dictionary. (disable class type checking so we can load
    // field)
    Info<< "Loading field " << fieldName << endl;
    const word oldTypeName = IOdictionary::typeName;
    const_cast<word&>(IOdictionary::typeName) = word::null;

    IOdictionary fieldDict
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    const_cast<word&>(IOdictionary::typeName) = oldTypeName;
    // Fake type back to what was in field
    const_cast<word&>(fieldDict.type()) = fieldDict.headerClassName();



    dictionary& boundaryField = fieldDict.subDict("boundaryField");

    label nChanged = 0;

    forAllConstIter(HashTable<word>, thisNames, iter)
    {
        const word& patchName = iter.key();
        const word& newName = iter();

        Info<< "Looking for entry for patch " << patchName << endl;

        // Find old patch name either direct or through wildcards
        // Find new patch name direct only

        if
        (
            boundaryField.found(patchName)
        && !boundaryField.found(newName, false, false)
        )
        {
            Info<< "    Changing entry " << patchName << " to " << newName
                << endl;

            dictionary& patchDict = boundaryField.subDict(patchName);

            if (patchDict.found("value"))
            {
                // Remove any value field since wrong size.
                patchDict.remove("value");
            }


            boundaryField.changeKeyword(patchName, newName);
            boundaryField.add
            (
                nbrNames[patchName],
                patchDict
            );
            Info<< "    Adding entry " << nbrNames[patchName] << endl;

            nChanged++;
        }
    }

    // Info<< "New boundaryField:" << boundaryField << endl;

    if (returnReduce(nChanged, sumOp<label>()) > 0)
    {
        if (isTestRun)
        {
            // Info<< "-test option: no changes made" << endl;
        }
        else
        {
            if (mvBak(fieldDict.objectPath(), "old"))
            {
                Info<< "Backup to    "
                    << (fieldDict.objectPath() + ".old") << nl;
            }

            Info<< "Write  to    "
                << fieldDict.objectPath() << endl;
            fieldDict.regIOobject::write();
        }
    }
    else
    {
        Info<< "No changes made to field " << fieldName << endl;
    }
    Info<< endl;
}


void rewriteFields
(
    const bool isTestRun,
    const Time& runTime,
    const wordList& fieldNames,
    const HashTable<word>& thisNames,
    const HashTable<word>& nbrNames
)
{
    forAll(fieldNames, i)
    {
        rewriteField
        (
            isTestRun,
            runTime,
            fieldNames[i],
            thisNames,
            nbrNames
        );
    }
}



int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    argList::addBoolOption("test", "test only; do not change any files");
    argList::addBoolOption
    (
        "enableFunctionEntries",
        "enable expansion of dictionary directives - #include, #codeStream etc"
    );
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"


    // Make sure we do not use the master-only reading since we read
    // fields (different per processor) as dictionaries.
    regIOobject::fileModificationChecking = regIOobject::timeStamp;


    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool isTestRun = args.optionFound("test");
    if (isTestRun)
    {
        Info<< "-test option: no changes made" << nl << endl;
    }
    const bool enableEntries = args.optionFound("enableFunctionEntries");


    Foam::word regionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", regionName);

    fileName regionPrefix = "";
    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }


    // Per cyclic patch the new name for this side and the other side
    HashTable<word> thisNames;
    HashTable<word> nbrNames;

    // Rewrite constant boundary file. Return any patches that have been split.
    IOobject io
    (
        "boundary",
        runTime.constant(),
        polyMesh::meshSubDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (io.typeHeaderOk<IOPtrList<entry>>(false))
    {
        rewriteBoundary
        (
            isTestRun,
            io,
            regionPrefix,
            thisNames,
            nbrNames
        );
    }



    // Convert any fields

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time: " << runTime.timeName() << endl;

        // See if mesh in time directory
        IOobject io
        (
            "boundary",
            runTime.timeName(),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (io.typeHeaderOk<IOPtrList<entry>>(false))
        {
            rewriteBoundary
            (
                isTestRun,
                io,
                regionPrefix,
                thisNames,
                nbrNames
            );
        }


        IOobjectList objects(runTime, runTime.timeName());


        int oldFlag = entry::disableFunctionEntries;
        if (!enableEntries)
        {
            // By default disable dictionary expansion for fields
            entry::disableFunctionEntries = 1;
        }

        // volFields
        // ~~~~~~~~~

        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(volScalarField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(volVectorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(volSphericalTensorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(volSymmTensorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(volTensorField::typeName),
            thisNames,
            nbrNames
        );


        // pointFields
        // ~~~~~~~~~~~

        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(pointScalarField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(pointVectorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(pointSphericalTensorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(pointSymmTensorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(pointTensorField::typeName),
            thisNames,
            nbrNames
        );


        // surfaceFields
        // ~~~~~~~~~~~

        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(surfaceScalarField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(surfaceVectorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(surfaceSphericalTensorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(surfaceSymmTensorField::typeName),
            thisNames,
            nbrNames
        );
        rewriteFields
        (
            isTestRun,
            runTime,
            objects.names(surfaceTensorField::typeName),
            thisNames,
            nbrNames
        );

        entry::disableFunctionEntries = oldFlag;
    }

    return 0;
}


// ************************************************************************* //
