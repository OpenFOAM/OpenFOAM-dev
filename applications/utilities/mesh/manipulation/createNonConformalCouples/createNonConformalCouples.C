/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    createNonConformalCouples

Description
    Utility to create non-conformal couples between non-coupled patches.

Usage
    \b createNonConformalCouples <patch1> <patch2>

    Options:
      - \par -overwrite \n
        Replace the old mesh with the new one, rather than writing the new one
        into a separate time directory

      - \par -fields \n
        Add non-conformal boundary conditions to the fields.

Note
    If run with two arguments, these arguments specify the patches between
    which a single couple is to be created. The resulting couple will not have
    a transformation.

Usage
    \b createNonConformalCouples

Note
    If run without arguments then settings are read from a \b
    system/createNonConformalCouplesDict dictionary (or from a different
    dictionary specified by the \b -dict option). This dictionary can specify
    the creation of multiple couples and/or couples with transformations.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMeshStitchersStationary.H"
#include "fvMeshTools.H"
#include "IOobjectList.H"
#include "nonConformalCyclicPolyPatch.H"
#include "nonConformalErrorPolyPatch.H"
#include "nonConformalProcessorCyclicPolyPatch.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "systemDict.H"
#include "Time.H"

#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    #include "addDictOption.H"

    const bool haveArgs = argList::hasArgs(argc, argv);
    if (haveArgs)
    {
        argList::validArgs.append("patch1");
        argList::validArgs.append("patch2");
        argList::addBoolOption
        (
            "fields",
            "add non-conformal boundary conditions to the fields"
        );
    }

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"

    // Flag to determine whether or not patches are added to fields
    bool fields;

    // Patch names between which to create couples, field dictionaries, the
    // associated cyclic name prefix and transformation (if any)
    List<Pair<word>> patchNames;
    List<Pair<dictionary>> patchFieldDicts;
    wordList cyclicNames;
    List<cyclicTransform> transforms;

    // If there are patch name arguments, then we assume fields are not being
    // changed, the cyclic name is just the cyclic typename, and that there is
    // no transformation. If there are no arguments then get all this
    // information from the system dictionary.
    if (haveArgs)
    {
        fields = args.optionFound("fields");

        patchNames.append(Pair<word>(args[1], args[2]));
        patchFieldDicts.append(Pair<dictionary>());
        cyclicNames.append(nonConformalCyclicPolyPatch::typeName);
        transforms.append(cyclicTransform(true));
    }
    else
    {
        static const word dictName("createNonConformalCouplesDict");

        IOdictionary dict(systemDict(dictName, args, runTime));

        fields = dict.lookupOrDefault<bool>("fields", false);

        const dictionary& couplesDict =
            dict.optionalSubDict("nonConformalCouples");

        forAllConstIter(dictionary, couplesDict, iter)
        {
            if (!iter().isDict()) continue;

            const dictionary& subDict = iter().dict();

            const bool havePatches = subDict.found("patches");
            const bool haveOwnerNeighbour =
                subDict.found("owner") || subDict.found("neighbour");

            if (havePatches == haveOwnerNeighbour)
            {
                FatalIOErrorInFunction(subDict)
                    << "Patches should be specified with either a single "
                    << "\"patches\" entry with a pair of patch names, or with "
                    << "two sub-dictionaries named \"owner\" and "
                    << "\"neighbour\"." << exit(FatalIOError);
            }

            if (havePatches)
            {
                patchNames.append(subDict.lookup<Pair<word>>("patches"));
                patchFieldDicts.append(Pair<dictionary>());
            }

            if (haveOwnerNeighbour)
            {
                const dictionary& ownerDict = subDict.subDict("owner");
                const dictionary& neighbourDict = subDict.subDict("neighbour");

                patchNames.append
                (
                    Pair<word>
                    (
                        ownerDict["patch"],
                        neighbourDict["patch"]
                    )
                );
                patchFieldDicts.append
                (
                    Pair<dictionary>
                    (
                        ownerDict.subOrEmptyDict("patchFields"),
                        neighbourDict.subOrEmptyDict("patchFields")
                    )
                );
            }

            cyclicNames.append(subDict.dictName());
            transforms.append(cyclicTransform(subDict, true));
        }
    }

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    #include "createNamedMesh.H"

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const bool overwrite = args.optionFound("overwrite");

    const word oldInstance = mesh.pointsInstance();

    // Read the fields
    IOobjectList objects(mesh, runTime.name());
    if (fields) Info<< "Reading geometric fields" << nl << endl;
    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"
    if (fields) Info<< endl;

    // Make sure the mesh is not connected before couples are added
    mesh.conform();

    // Start building lists of patches and patch-fields to add
    List<polyPatch*> newPatches;
    List<dictionary> newPatchFieldDicts;

    // Find the first processor patch and face
    label firstProcPatchi = patches.size(), firstProcFacei = mesh.nFaces();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<processorPolyPatch>(pp) && firstProcPatchi == patches.size())
        {
            firstProcPatchi = patchi;
            firstProcFacei = pp.start();
        }

        if (!isA<processorPolyPatch>(pp) && firstProcPatchi != patches.size())
        {
            FatalErrorInFunction
                << "Processor patches do not follow boundary patches"
                << exit(FatalError);
        }
    }

    // Clone the non-processor patches
    for (label patchi = 0; patchi < firstProcPatchi; ++ patchi)
    {
        const polyPatch& pp = patches[patchi];

        newPatches.append
        (
            pp.clone(patches, patchi, pp.size(), pp.start()).ptr()
        );
        newPatchFieldDicts.append
        (
            dictionary()
        );
    }

    // Convenience function to generate patch names for the owner or neighbour
    auto nccPatchNames = [&](const label i)
    {
        return
            Pair<word>
            (
                cyclicNames[i] + "_on_" + patchNames[i][0],
                cyclicNames[i] + "_on_" + patchNames[i][1]
            );
    };

    // Add the cyclic patches
    forAll(patchNames, i)
    {
        Info<< indent << "Adding "
            << nonConformalCyclicPolyPatch::typeName
            << " interfaces between patches: " << incrIndent << nl
            << indent << patchNames[i] << decrIndent << nl
            << indent << "Named:" << incrIndent << nl
            << indent << Pair<word>(nccPatchNames(i)) << decrIndent << nl
            << indent << "With transform: " << incrIndent << nl;
        transforms[i].write(Info);
        Info<< decrIndent << nl;

        newPatches.append
        (
            new nonConformalCyclicPolyPatch
            (
                nccPatchNames(i)[0],
                0,
                firstProcFacei,
                newPatches.size(),
                patches,
                nonConformalCyclicPolyPatch::typeName,
                nccPatchNames(i)[1],
                patchNames[i][0],
                transforms[i]
            )
        );
        newPatchFieldDicts.append
        (
            patchFieldDicts[i][0]
        );
        newPatches.append
        (
            new nonConformalCyclicPolyPatch
            (
                nccPatchNames(i)[1],
                0,
                firstProcFacei,
                newPatches.size(),
                patches,
                nonConformalCyclicPolyPatch::typeName,
                nccPatchNames(i)[0],
                patchNames[i][1],
                inv(transforms[i])
            )
        );
        newPatchFieldDicts.append
        (
            patchFieldDicts[i][1]
        );
    }

    // Add the error patches. Note there is only one for each source patch,
    // regardless of how many interfaces are attached to that patch.
    auto appendErrorPatches = [&](const bool owner)
    {
        wordHashSet patchANames;
        forAll(patchNames, i)
        {
            patchANames.insert(patchNames[i][!owner]);
        }
        forAllConstIter(wordHashSet, patchANames, iter)
        {
            newPatches.append
            (
                new nonConformalErrorPolyPatch
                (
                    nonConformalErrorPolyPatch::typeName + "_on_" + iter.key(),
                    0,
                    firstProcFacei,
                    newPatches.size(),
                    patches,
                    nonConformalErrorPolyPatch::typeName,
                    iter.key()
                )
            );
            newPatchFieldDicts.append
            (
                dictionary()
            );
        }
    };
    appendErrorPatches(true);
    appendErrorPatches(false);

    // Clone the processor patches
    for (label patchi = firstProcPatchi; patchi < patches.size(); ++ patchi)
    {
        const polyPatch& pp = patches[patchi];

        newPatches.append
        (
            pp.clone(patches, newPatches.size(), pp.size(), pp.start()).ptr()
        );
        newPatchFieldDicts.append
        (
            dictionary()
        );
    }

    // Add the processor cyclic patches
    if (Pstream::parRun())
    {
        forAll(patchNames, i)
        {
            const polyPatch& patch1 = patches[patchNames[i][0]];
            const polyPatch& patch2 = patches[patchNames[i][1]];

            boolList procHasPatch1(Pstream::nProcs(), false);
            procHasPatch1[Pstream::myProcNo()] = !patch1.empty();
            Pstream::gatherList(procHasPatch1);
            Pstream::scatterList(procHasPatch1);

            boolList procHasPatch2(Pstream::nProcs(), false);
            procHasPatch2[Pstream::myProcNo()] = !patch2.empty();
            Pstream::gatherList(procHasPatch2);
            Pstream::scatterList(procHasPatch2);

            // Multiple cyclic interfaces must be ordered in a specific way for
            // processor communication to function correctly.
            //
            // A communication that is sent from the cyclic owner is received
            // on the cyclic neighbour and vice versa. Therefore, in a coupled
            // pair of processors if one sends the owner first the other must
            // receive the neighbour first.
            //
            // We ensure the above by ordering the patches so that for the
            // lower indexed processor the owner interface comes first, and for
            // the higher indexed processor the neighbour comes first.

            auto appendProcPatches = [&](const bool owner, const bool first)
            {
                const boolList& procHasPatchA =
                    owner ? procHasPatch1 : procHasPatch2;
                const boolList& procHasPatchB =
                    owner ? procHasPatch2 : procHasPatch1;

                if (procHasPatchA[Pstream::myProcNo()])
                {
                    forAll(procHasPatchB, proci)
                    {
                        if
                        (
                            (
                                (first && proci < Pstream::myProcNo())
                             || (!first && proci > Pstream::myProcNo())
                            )
                         && procHasPatchB[proci]
                        )
                        {
                            newPatches.append
                            (
                                new nonConformalProcessorCyclicPolyPatch
                                (
                                    0,
                                    mesh.nFaces(),
                                    newPatches.size(),
                                    patches,
                                    Pstream::myProcNo(),
                                    proci,
                                    nccPatchNames(i)[!owner],
                                    patchNames[i][!owner]
                                )
                            );
                            newPatchFieldDicts.append
                            (
                                patchFieldDicts[i][!owner]
                            );
                        }
                    }
                }
            };

            appendProcPatches(true, true);
            appendProcPatches(false, true);
            appendProcPatches(false, false);
            appendProcPatches(true, false);
        }
    }

    // Re-patch the mesh. Note that new patches are all constraints, so the
    // dictionary and patch type do not get used. Overrides will be handled
    // later, once all patches have been added and the mesh has been stitched.
    forAll(newPatches, newPatchi)
    {
        fvMeshTools::addPatch
        (
            mesh,
            *newPatches[newPatchi],
            dictionary(),
            calculatedFvPatchField<scalar>::typeName,
            false
        );
    }

    // Connect the mesh so that the new stitching topology gets written out
    fvMeshStitchers::stationary(mesh).connect(false, false, false);

    // Set the fields on the new patches. All new patches are constraints, so
    // this should only be creating overrides; e.g., jump cyclics.
    forAll(newPatches, newPatchi)
    {
        if (!newPatchFieldDicts[newPatchi].empty())
        {
            fvMeshTools::setPatchFields
            (
                mesh,
                newPatchi,
                newPatchFieldDicts[newPatchi]
            );
        }
    }

    mesh.setInstance(runTime.name());

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info<< "Writing mesh to " << runTime.name() << nl << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
