/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    changeDictionary

Description
    Utility to change dictionary entries, e.g. can be used to change the patch
    type in the field and polyMesh/boundary files.

    Reads dictionaries (fields) and entries to change from a dictionary.
    E.g. to make the \em movingWall a \em fixedValue for \em p but all other
    \em Walls a zeroGradient boundary condition, the
    \c system/changeDictionaryDict would contain the following:
    \verbatim
    p                           // field to change
    {
        boundaryField
        {
            ".*Wall"            // entry to change
            {
                type            zeroGradient;
            }
            movingWall          // entry to change
            {
                type            fixedValue;
                value           uniform 123.45;
            }
        }
    }
    \endverbatim
    Replacement entries starting with '~' will remove the entry.

Usage
    \b changeDictionary [OPTION]

    Options:
      - \par -subDict
        Specify the subDict name of the replacements dictionary.

      - \par -literalRE
        Do not interpret regular expressions; treat them as any other keyword.

      - \par -enableFunctionEntries
        Enable function entries (default: disabled)

      - \par -disablePatchGroups
        Disable the default checking for keys being patchGroups

    Note:
        changeDictionary has been superseded by foamDictionary
        and is now deprecated.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "volFields.H"
#include "stringListOps.H"
#include "timeSelector.H"
#include "systemDict.H"
#include "mergeDictionaries.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Extract groupPatch info from boundary file info
HashTable<wordList, word> extractPatchGroups(const dictionary& boundaryDict)
{
    HashTable<wordList, word> groupToPatch;

    forAllConstIter(dictionary, boundaryDict, iter)
    {
        const word& patchName = iter().keyword();
        const dictionary& patchDict = iter().dict();

        wordList groups;
        if (patchDict.readIfPresent("inGroups", groups))
        {
            forAll(groups, i)
            {
                HashTable<wordList, word>::iterator fndGroup = groupToPatch.find
                (
                    groups[i]
                );
                if (fndGroup == groupToPatch.end())
                {
                    groupToPatch.insert(groups[i], wordList(1, patchName));
                }
                else
                {
                    fndGroup().append(patchName);
                }
            }
        }
    }
    return groupToPatch;
}


int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    argList::addOption
    (
        "subDict",
        "name",
        "specify the subDict name of the replacements dictionary"
    );
    argList::addOption
    (
        "instance",
        "name",
        "override instance setting (default is the time name)"
    );

    // Add explicit time option
    timeSelector::addOptions();

    argList::addBoolOption
    (
        "literalRE",
        "treat regular expressions literally (i.e., as a keyword)"
    );
    argList::addBoolOption
    (
        "enableFunctionEntries",
        "enable expansion of dictionary directives - #include, #codeStream etc"
    );
    argList::addBoolOption
    (
        "disablePatchGroups",
        "disable matching keys to patch groups"
    );

    #include "addMeshOption.H"
    #include "addRegionOption.H"

    #include "setRootCase.H"

    Warning
        << nl
        << "changeDictionary has been superseded by foamDictionary"
           " and is now deprecated."
        << nl << endl;

    #include "createTime.H"

    // Optionally override controlDict time with -time options
    instantList times = timeSelector::selectIfPresent(runTime, args);
    if (times.size() < 1)
    {
        FatalErrorInFunction
            << "No times selected." << exit(FatalError);
    }
    forAll(times, timei)
    {
        word instance;
        if (args.optionFound("instance"))
        {
            if (times.size() > 1)
            {
                FatalErrorInFunction
                    << "Multiple times selected with 'instance' option"
                    << exit(FatalError);
            }

            args.optionLookup("instance")() >> instance;
        }
        else
        {
            runTime.setTime(times[timei], timei);
            instance = runTime.name();
        }

        #include "createSpecifiedMeshNoChangers.H"

        const bool literalRE = args.optionFound("literalRE");
        if (literalRE)
        {
            Info<< "Not interpreting any regular expressions (RE)"
                << " in the changeDictionaryDict." << endl
                << "Instead they are handled as any other entry, i.e. added if"
                << " not present." << endl;
        }

        const bool enableEntries = args.optionFound("enableFunctionEntries");
        if (enableEntries)
        {
            Info<< "Allowing dictionary preprocessing "
                   "('#include', '#codeStream')."
                << endl;
        }

        int oldFlag = entry::disableFunctionEntries;
        if (!enableEntries)
        {
            // By default disable dictionary expansion for fields
            entry::disableFunctionEntries = 1;
        }


        const bool disablePatchGroups = args.optionFound("disablePatchGroups");
        if (disablePatchGroups)
        {
            Info<< "Not interpreting any keys in the changeDictionary"
                << " as patchGroups"
                << endl;
        }


        fileName regionPrefix = "";
        if (regionName != fvMesh::defaultRegion)
        {
            regionPrefix = regionName;
        }


        // Make sure we do not use the master-only reading since we read
        // fields (different per processor) as dictionaries.
        regIOobject::fileModificationChecking = regIOobject::timeStamp;


        // Get the replacement rules from a dictionary
        const dictionary dict(systemDict("changeDictionaryDict", args, mesh));

        const dictionary* replaceDictsPtr = &dict;

        if (args.optionFound("subDict"))
        {
            word subDictName(args.optionLookup("subDict")());
            replaceDictsPtr = &dict.subDict(subDictName);
        }

        const dictionary& replaceDicts = *replaceDictsPtr;

        Info<< "Read dictionary " << dict.name()
            << " with replacements for dictionaries "
            << replaceDicts.toc() << endl;



        // Always read boundary to get patch groups
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Reading polyMesh/boundary file to extract patch names"
            << endl;

        // Read PtrList of dictionary as dictionary.
        const word oldTypeName = IOPtrList<entry>::typeName;
        const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
        IOPtrList<entry> dictList
        (
            IOobject
            (
                "boundary",
                runTime.findInstance
                (
                    regionPrefix/polyMesh::meshSubDir,
                    "boundary",
                    IOobject::READ_IF_PRESENT
                ),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );
        const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;

        // Fake type back to what was in field
        const_cast<word&>(dictList.type()) = dictList.headerClassName();

        // Temporary convert to dictionary
        dictionary fieldDict;
        forAll(dictList, i)
        {
            fieldDict.add(dictList[i].keyword(), dictList[i].dict());
        }

        if (dictList.size())
        {
            Info<< "Loaded dictionary " << dictList.name()
                << " with entries " << fieldDict.toc() << endl;
        }

        // Extract any patchGroups information (= shortcut for set of
        // patches)
        HashTable<wordList, word> patchGroups;
        if (!disablePatchGroups)
        {
            patchGroups = extractPatchGroups(fieldDict);
            if (patchGroups.size())
            {
                Info<< "Extracted patch groups:" << endl;
                wordList groups(patchGroups.sortedToc());
                forAll(groups, i)
                {
                    Info<< "    group " << groups[i] << " with patches "
                        << patchGroups[groups[i]] << endl;
                }
            }
        }


        // Every replacement is a dictionary name and a keyword in this

        forAllConstIter(dictionary, replaceDicts, fieldIter)
        {
            const word& fieldName = fieldIter().keyword();
            Info<< "Replacing entries in dictionary " << fieldName << endl;

            // Handle 'boundary' specially:
            // - is PtrList of dictionaries
            // - is in polyMesh/
            if (fieldName == "boundary")
            {
                Info<< "Special handling of " << fieldName
                    << " as polyMesh/boundary file." << endl;

                // Get the replacement dictionary for the field
                const dictionary& replaceDict = fieldIter().dict();
                Info<< "Merging entries from " << replaceDict.toc() << endl;

                // Merge the replacements in
                mergeDictionaries
                (
                    fieldDict,
                    replaceDict,
                    !literalRE,
                    patchGroups
                );

                Info<< "fieldDict:" << fieldDict << endl;

                // Convert back into dictList
                wordList doneKeys(dictList.size());

                label nEntries = fieldDict.size();

                forAll(dictList, i)
                {
                    doneKeys[i] = dictList[i].keyword();
                    dictList.set
                    (
                        i,
                        fieldDict.lookupEntry
                        (
                            doneKeys[i],
                            false,
                            true
                        ).clone()
                    );
                    fieldDict.remove(doneKeys[i]);
                }

                // Add remaining entries
                label sz = dictList.size();
                dictList.setSize(nEntries);
                forAllConstIter(dictionary, fieldDict, iter)
                {
                    dictList.set(sz++, iter().clone());
                }

                Info<< "Writing modified " << fieldName << endl;
                dictList.writeObject
                (
                    runTime.writeFormat(),
                    runTime.writeFormat(),
                    IOstream::UNCOMPRESSED,
                    true
                );
            }
            else
            {
                // Read dictionary. (disable class type checking so we can load
                // field)
                Info<< "Loading dictionary " << fieldName << endl;
                const word oldTypeName = localIOdictionary::typeName;
                const_cast<word&>(localIOdictionary::typeName) = word::null;

                localIOdictionary fieldDict
                (
                    IOobject
                    (
                        fieldName,
                        instance,
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE,
                        false
                    )
                );

                const_cast<word&>(localIOdictionary::typeName) = oldTypeName;

                // Fake type back to what was in field
                const_cast<word&>(fieldDict.type()) =
                    fieldDict.headerClassName();

                Info<< "Loaded dictionary " << fieldName
                    << " with entries " << fieldDict.toc() << endl;

                // Get the replacement dictionary for the field
                const dictionary& replaceDict = fieldIter().dict();
                Info<< "Merging entries from " << replaceDict.toc() << endl;

                // Merge the replacements in
                mergeDictionaries
                (
                    fieldDict,
                    replaceDict,
                    !literalRE,
                    patchGroups
                );

                Info<< "Writing modified fieldDict " << fieldName << endl;
                fieldDict.regIOobject::write();
            }
        }

        entry::disableFunctionEntries = oldFlag;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
