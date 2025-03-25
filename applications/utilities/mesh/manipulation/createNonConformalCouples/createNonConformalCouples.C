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
    createNonConformalCouples

Description
    Utility to create non-conformal couples between non-coupled patches.

Usage
    \b createNonConformalCouples \<patch1\> \<patch2\>

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
#include "fvMeshTools.H"
#include "hashedWordList.H"
#include "IOobjectList.H"
#include "MultiRegionList.H"
#include "nonConformalCyclicPolyPatch.H"
#include "nonConformalErrorPolyPatch.H"
#include "nonConformalMappedWallPolyPatch.H"
#include "nonConformalProcessorCyclicPolyPatch.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "stationary_fvMeshStitcher.H"
#include "systemDict.H"
#include "Time.H"

#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct nonConformalCouple
{
    // Public Data

        Pair<word> regionNames;

        Pair<word> origPatchNames;

        word ncPatchType;

        Pair<word> ncPatchNames;

        Pair<dictionary> ncPatchDicts;

        Pair<dictionary> ncPatchFieldDicts;

        cyclicTransform transform;


    // Constructors

        // Construct null (for List)
        nonConformalCouple()
        {}

        // Construct from arguments
        nonConformalCouple
        (
            const word& primaryRegionName,
            const argList& args
        )
        :
            regionNames(primaryRegionName, primaryRegionName),
            origPatchNames(args[1], args[2]),
            ncPatchType(nonConformalCyclicPolyPatch::typeName),
            ncPatchNames
            (
                ncPatchType + "_on_" + args[1],
                ncPatchType + "_on_" + args[2]
            ),
            ncPatchDicts(),
            ncPatchFieldDicts(),
            transform(true)
        {}

        //- Construct from dictionary
        nonConformalCouple
        (
            const word& primaryRegionName,
            const dictionary& dict
        )
        {
            const bool haveRegion = dict.found("region");
            const bool haveRegions = dict.found("regions");

            if (haveRegion && haveRegions)
            {
                FatalIOErrorInFunction(dict)
                    << "Regions should be specified either by a \"region\" "
                    << "entry with a single region name (for in-region "
                    << "cyclics), or by a \"regions\" entry with a pair of "
                    << "names (for inter-region mapped walls)"
                    << exit(FatalIOError);
            }

            const bool havePatches = dict.found("patches");
            const bool haveOwnerNeighbour =
                dict.found("owner") || dict.found("neighbour");

            if (havePatches == haveOwnerNeighbour)
            {
                FatalIOErrorInFunction(dict)
                    << "Patches should be specified with either a single "
                    << "\"patches\" entry with a pair of patch names, "
                    << "or with two sub-dictionaries named \"owner\" and "
                    << "\"neighbour\"" << exit(FatalIOError);
            }

            if (havePatches)
            {
                const word thisRegionName =
                    haveRegion ? dict.lookup<word>("region") : word::null;

                regionNames =
                    haveRegion ? Pair<word>(thisRegionName, thisRegionName)
                  : haveRegions ? dict.lookup<Pair<word>>("regions")
                  : Pair<word>(primaryRegionName, primaryRegionName);
                origPatchNames =
                    dict.lookup<Pair<word>>("patches");
                ncPatchType =
                    dict.lookupOrDefault<word>
                    (
                        "type",
                        regionNames.first() == regionNames.second()
                      ? nonConformalCyclicPolyPatch::typeName
                      : nonConformalMappedWallPolyPatch::typeName
                    );
                ncPatchNames =
                    dict.lookupOrDefault<Pair<word>>
                    (
                        "names",
                        Pair<word>
                        (
                            dict.dictName() + "_on_" + origPatchNames[0],
                            dict.dictName() + "_on_" + origPatchNames[1]
                        )
                    );
                forAll(ncPatchDicts, i)
                {
                    ncPatchDicts[i] = dict;
                    ncPatchDicts[i].remove("region");
                    ncPatchDicts[i].remove("regions");
                    ncPatchDicts[i].remove("patches");
                    ncPatchDicts[i].remove("type");
                    ncPatchDicts[i].remove("names");
                    ncPatchDicts[i].remove(cyclicTransform::keywords);
                }
                ncPatchFieldDicts = Pair<dictionary>();
                transform = cyclicTransform(dict, true);
            }
            else
            {
                const dictionary& ownerDict = dict.subDict("owner");
                const dictionary& neighbourDict = dict.subDict("neighbour");

                regionNames =
                    Pair<word>
                    (
                        ownerDict.lookupOrDefault<word>
                        (
                            "region",
                            primaryRegionName
                        ),
                        neighbourDict.lookupOrDefault<word>
                        (
                            "region",
                            primaryRegionName
                        )
                    );
                origPatchNames =
                    Pair<word>
                    (
                        ownerDict.lookup<word>("patch"),
                        neighbourDict.lookup<word>("patch")
                    );
                ncPatchType =
                    dict.lookupOrDefault<word>
                    (
                        "type",
                        regionNames.first() == regionNames.second()
                      ? nonConformalCyclicPolyPatch::typeName
                      : nonConformalMappedWallPolyPatch::typeName
                    );
                ncPatchNames =
                    Pair<word>
                    (
                        ownerDict.lookupOrDefault<word>
                        (
                            "name",
                            dict.dictName() + "_on_" + origPatchNames[0]
                        ),
                        neighbourDict.lookupOrDefault<word>
                        (
                            "name",
                            dict.dictName() + "_on_" + origPatchNames[1]
                        )
                    );
                ncPatchDicts[0] = ownerDict;
                ncPatchDicts[1] = neighbourDict;
                forAll(ncPatchDicts, i)
                {
                    ncPatchDicts[i].remove("region");
                    ncPatchDicts[i].remove("patch");
                    ncPatchDicts[i].remove("name");
                    ncPatchDicts[i].remove("patchFields");
                }
                ncPatchFieldDicts =
                    Pair<dictionary>
                    (
                        ownerDict.subOrEmptyDict("patchFields"),
                        neighbourDict.subOrEmptyDict("patchFields")
                    );
                transform = cyclicTransform(dict, true);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void evaluateNonConformalProcessorCyclics(const fvMesh& mesh)
{
    UPtrList<VolField<Type>> fields(mesh.fields<VolField<Type>>());

    forAll(fields, i)
    {
        const label nReq = Pstream::nRequests();

        forAll(mesh.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                fields[i].boundaryFieldRef()[patchi];

            if (isA<nonConformalProcessorCyclicPolyPatch>(pf.patch().patch()))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mesh.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                fields[i].boundaryFieldRef()[patchi];

            if (isA<nonConformalProcessorCyclicPolyPatch>(pf.patch().patch()))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addMeshOption.H"
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
    #include "setMeshPath.H"
    #include "createTimeNoFunctionObjects.H"

    const Foam::word primaryRegionName =
        args.optionLookupOrDefault("region", polyMesh::defaultRegion);

    // Flag to determine whether or not patches are added to fields
    bool fields;

    // Read the couples to be created from arguments or a system dictionary
    List<nonConformalCouple> couples;
    if (haveArgs)
    {
        fields = args.optionFound("fields");

        couples.append(nonConformalCouple(primaryRegionName, args));
    }
    else
    {
        const dictionary dict
        (
            systemDict
            (
                "createNonConformalCouplesDict",
                args,
                runTime,
                primaryRegionName
            )
        );

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
                    << "\"patches\" entry with a pair of patch names, "
                    << "or with two sub-dictionaries named \"owner\" and "
                    << "\"neighbour\"." << exit(FatalIOError);
            }

            couples.append(nonConformalCouple(primaryRegionName, subDict));
        }
    }

    // Determine if operating on meshes within a sub-path
    const bool haveMeshPath = meshPath != word::null;

    // Create all the meshes needed
    hashedWordList regionNames;
    PtrList<fvMesh> regionMeshes;
    forAll(couples, i)
    {
        forAll(couples[i].regionNames, sidei)
        {
            const word& regionName = couples[i].regionNames[sidei];

            if (regionNames.found(regionName)) continue;

            const IOobject regionMeshIo
            (
                regionName,
                runTime.name(),
                meshPath,
                runTime,
                Foam::IOobject::MUST_READ
            );

            // If we are creating couples on a mesh within a mesh path
            // sub-directory then these couples will not be stitched so loading
            // the neighbouring regions is optional
            if (haveMeshPath && !polyMesh::found(regionMeshIo)) continue;

            regionNames.append(regionName);
            regionMeshes.append(new fvMesh(regionMeshIo, false));
        }
    }

    // Early exit if there is nothing to do
    if (regionMeshes.empty())
    {
        Info<< "Nothing to be done" << nl << endl
            << "End" << nl << endl;

        return 0;
    }

    const bool overwrite = args.optionFound("overwrite");

    const word oldInstance = regionMeshes[0].pointsInstance();

    // Read the fields
    if (fields) Info<< "Reading geometric fields" << nl << endl;

    #define DeclareRegionGeoTypeFields(Type, Geo)                              \
        PtrList<PtrList<Geo##Field<Type>>>                                     \
            CAT4(region, Geo, CAPITALIZE(Type), Fields)(regionMeshes.size());  \
        forAll(regionMeshes, regioni)                                          \
        {                                                                      \
            CAT4(region, Geo, CAPITALIZE(Type), Fields).set                    \
            (                                                                  \
                regioni,                                                       \
                new PtrList<Geo##Field<Type>>                                  \
            );                                                                 \
        }
    FOR_ALL_FIELD_TYPES(DeclareRegionGeoTypeFields, Vol);
    FOR_ALL_FIELD_TYPES(DeclareRegionGeoTypeFields, Surface);
    FOR_ALL_FIELD_TYPES(DeclareRegionGeoTypeFields, Point);
    #undef DeclareRegionGeoTypeFields

    if (fields)
    {
        MultiRegionUList<fvMesh> multiRegionMeshes(regionMeshes, false);

        forAll(regionMeshes, regioni)
        {
            RegionRef<fvMesh> mesh = multiRegionMeshes[regioni];

            IOobjectList objects(mesh, runTime.name());

            #define ReadRegionGeoTypeFields(Type, Geo, mesh)                   \
                ReadFields                                                     \
                (                                                              \
                    mesh,                                                      \
                    objects,                                                   \
                    CAT4(region, Geo, CAPITALIZE(Type), Fields)[regioni]       \
                );
            FOR_ALL_FIELD_TYPES(ReadRegionGeoTypeFields, Vol, mesh);
            FOR_ALL_FIELD_TYPES(ReadRegionGeoTypeFields, Surface, mesh);
            const pointMesh& pMesh = pointMesh::New(mesh);
            FOR_ALL_FIELD_TYPES(ReadRegionGeoTypeFields, Point, pMesh);
            #undef ReadRegionGeoTypeFields
        }

        Info<< endl;
    }

    if (!overwrite)
    {
        runTime++;
    }

    // Make sure the meshes are not connected before couples are added
    forAll(regionMeshes, regioni)
    {
        regionMeshes[regioni].conform();
    }

    // Find the first processor patch and face
    labelList regionFirstProcPatchis(regionMeshes.size());
    labelList regionFirstProcFaceis(regionMeshes.size());;
    forAll(regionMeshes, regioni)
    {
        const fvMesh& mesh = regionMeshes[regioni];
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        label& firstProcPatchi = regionFirstProcPatchis[regioni];
        label& firstProcFacei = regionFirstProcFaceis[regioni];

        firstProcPatchi = patches.size();
        firstProcFacei = mesh.nFaces();

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            const bool isProcPp = isA<processorPolyPatch>(pp);

            if (isProcPp && firstProcPatchi == patches.size())
            {
                firstProcPatchi = patchi;
                firstProcFacei = pp.start();
            }

            if (!isProcPp && firstProcPatchi != patches.size())
            {
                FatalErrorInFunction
                    << "Processor patches of region " << regionNames[regioni]
                    << " do not follow boundary patches"
                    << exit(FatalError);
            }
        }
    }

    // Start building lists of patches and patch-fields to add
    List<List<polyPatch*>> newPatches(regionMeshes.size());
    List<boolList> newPatchIsCouple(regionMeshes.size());
    List<List<dictionary>> newPatchFieldDicts(regionMeshes.size());

    // Clone the non-processor patches
    forAll(regionMeshes, regioni)
    {
        const fvMesh& mesh = regionMeshes[regioni];
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const label firstProcPatchi = regionFirstProcPatchis[regioni];

        for (label patchi = 0; patchi < firstProcPatchi; ++ patchi)
        {
            const polyPatch& pp = patches[patchi];

            newPatches[regioni].append
            (
                pp.clone
                (
                    patches,
                    patchi,
                    pp.size(),
                    pp.start()
                ).ptr()
            );
            newPatchIsCouple[regioni].append(false);
            newPatchFieldDicts[regioni].append(dictionary());
        }
    }

    // Add the non-processor coupled patches
    forAll(couples, couplei)
    {
        const nonConformalCouple& couple = couples[couplei];

        Info<< indent << "Adding " << couple.ncPatchType
            << " interfaces between patches: " << incrIndent << nl
            << indent << couple.origPatchNames << decrIndent << nl
            << "In regions: " << incrIndent << nl
            << indent << couple.regionNames << decrIndent << nl
            << indent << "Named:" << incrIndent << nl
            << indent << couple.ncPatchNames << decrIndent << nl
            << indent << "With transform: " << incrIndent << nl;
        couple.transform.write(Info);
        Info<< decrIndent << nl;

        auto appendPatch = [&](const bool owner)
        {
            if (!regionNames.found(couple.regionNames[!owner])) return;

            const label regioni =
                regionNames[couple.regionNames[!owner]];

            dictionary patchDict
            (
                "type", couple.ncPatchType,
                "nFaces", 0,
                "startFace", regionFirstProcFaceis[regioni],
                "originalPatch", couple.origPatchNames[!owner],
                "neighbourRegion", couple.regionNames[owner],
                "neighbourPatch", couple.ncPatchNames[owner],
                "owner", owner
            );

            {
                OStringStream oss;
                (owner ? couple.transform : inv(couple.transform)).write(oss);
                patchDict.merge(IStringStream(oss.str())());
            }

            patchDict.merge(couple.ncPatchDicts[!owner]);

            newPatches[regioni].append
            (
                polyPatch::New
                (
                    couple.ncPatchNames[!owner],
                    patchDict,
                    newPatches[regioni].size(),
                    regionMeshes[regioni].boundaryMesh()
                ).ptr()
            );
            newPatchIsCouple[regioni].append(true);
            newPatchFieldDicts[regioni].append
            (
                couple.ncPatchFieldDicts[!owner]
            );
        };
        appendPatch(true);
        appendPatch(false);
    }

    // Add the error patches. Note there is only one for each original patch,
    // regardless of how many couplings are attached to that patch.
    {
        // Create a table of unique original patches
        HashTable<label, Pair<word>, Hash<Pair<word>>> regionOrigPatchToIndex;
        DynamicList<Pair<word>> regionOrigPatches;
        forAll(couples, couplei)
        {
            const nonConformalCouple& couple = couples[couplei];
            forAll(couple.regionNames, i)
            {
                const Pair<word> regionOrigPatch
                (
                    couple.regionNames[i],
                    couple.origPatchNames[i]
                );

                if (!regionOrigPatchToIndex.found(regionOrigPatch))
                {
                    regionOrigPatchToIndex.insert
                    (
                        regionOrigPatch,
                        regionOrigPatches.size()
                    );
                    regionOrigPatches.append(regionOrigPatch);
                }
            }
        }

        // Add an error patch for each unique original patch
        forAll(regionOrigPatches, i)
        {
            const word& regionName = regionOrigPatches[i].first();

            if (!regionNames.found(regionName)) continue;

            const label regioni = regionNames[regionName];
            const word& origPatchName = regionOrigPatches[i].second();

            newPatches[regioni].append
            (
                new nonConformalErrorPolyPatch
                (
                    nonConformalErrorPolyPatch::typeName
                  + "_on_"
                  + origPatchName,
                    0,
                    regionFirstProcFaceis[regioni],
                    newPatches[regioni].size(),
                    regionMeshes[regioni].boundaryMesh(),
                    nonConformalErrorPolyPatch::typeName,
                    origPatchName
                )
            );
            newPatchIsCouple[regioni].append(false);
            newPatchFieldDicts[regioni].append(dictionary());
        }
    }

    // Clone the processor patches
    forAll(regionMeshes, regioni)
    {
        const fvMesh& mesh = regionMeshes[regioni];
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const label firstProcPatchi = regionFirstProcPatchis[regioni];

        for (label patchi = firstProcPatchi; patchi < patches.size(); ++ patchi)
        {
            const polyPatch& pp = patches[patchi];

            newPatches[regioni].append
            (
                pp.clone
                (
                    patches,
                    newPatches[regioni].size(),
                    pp.size(),
                    pp.start()
                ).ptr()
            );
            newPatchIsCouple[regioni].append(false);
            newPatchFieldDicts[regioni].append(dictionary());
        }
    }

    // Add the processor cyclic patches
    if (Pstream::parRun())
    {
        forAll(couples, couplei)
        {
            const nonConformalCouple& couple = couples[couplei];

            if (couple.ncPatchType != nonConformalCyclicPolyPatch::typeName)
                continue;

            if (!regionNames.found(couples[couplei].regionNames[0]))
                continue;

            const label regioni = regionNames[couples[couplei].regionNames[0]];
            const fvMesh& mesh = regionMeshes[regioni];
            const polyBoundaryMesh& patches = mesh.boundaryMesh();

            const polyPatch& patch1 = patches[couple.origPatchNames[0]];
            const polyPatch& patch2 = patches[couple.origPatchNames[1]];

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
                                (first && proci > Pstream::myProcNo())
                             || (!first && proci < Pstream::myProcNo())
                            )
                         && procHasPatchB[proci]
                        )
                        {
                            newPatches[regioni].append
                            (
                                new nonConformalProcessorCyclicPolyPatch
                                (
                                    0,
                                    mesh.nFaces(),
                                    newPatches.size(),
                                    patches,
                                    Pstream::myProcNo(),
                                    proci,
                                    couple.ncPatchNames[!owner],
                                    couple.origPatchNames[!owner]
                                )
                            );
                            newPatchIsCouple[regioni].append(true);
                            newPatchFieldDicts[regioni].append
                            (
                                couple.ncPatchFieldDicts[!owner]
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

    // Re-patch the meshes. Create constraint or calculated patch fields at
    // this stage; don't apply the patch field dictionaries. The patch fields
    // will be specified properly later after all patches have been added and
    // the meshes have been stitched. That way the geometry is fully available
    // for construction of the patch fields.
    forAll(regionMeshes, regioni)
    {
        forAll(newPatches[regioni], newPatchi)
        {
            fvMeshTools::addPatch
            (
                regionMeshes[regioni],
                *newPatches[regioni][newPatchi]
            );
        }
    }

    // Connect the meshes
    {
        MultiRegionUList<fvMesh> multiRegionMeshes(regionMeshes, false);
        forAll(regionMeshes, regioni)
        {
            RegionRef<fvMesh> mesh = multiRegionMeshes[regioni];
            fvMeshStitchers::stationary(mesh).connect(false, false, false);
        }

        if (!haveMeshPath) Info<< endl;
    }

    // Set the fields on the new patches. This allows constraint types (e.g.,
    // nonConformalCyclic) to be set by default, but requires a dictionary
    // setting for non-constrained types (e.g., nonConformalMappedWall). It
    // also allows constraint types to be overridden (e.g., with jumpCyclic) if
    // there is a dictionary present.
    forAll(regionMeshes, regioni)
    {
        forAll(newPatches[regioni], newPatchi)
        {
            if (newPatchIsCouple[regioni][newPatchi])
            {
                fvMeshTools::setPatchFields
                (
                    regionMeshes[regioni],
                    newPatchi,
                    newPatchFieldDicts[regioni][newPatchi]
                );
            }
        }
    }

    // Communicate values across non-conformal processor cyclics so that they
    // contain valid values that can be written to disk
    if (Pstream::parRun())
    {
        forAll(regionMeshes, regioni)
        {
            const fvMesh& mesh = regionMeshes[regioni];

            #define EVALUATE_NON_CONFORMAL_PROCESSOR_CYCLICS(Type, nullArg) \
                evaluateNonConformalProcessorCyclics<Type>(mesh);
            FOR_ALL_FIELD_TYPES(EVALUATE_NON_CONFORMAL_PROCESSOR_CYCLICS)
            #undef EVALUATE_NON_CONFORMAL_PROCESSOR_CYCLICS
        }
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    // Set the instance so that the mesh writes
    const word newInstance = overwrite ? oldInstance : runTime.name();
    forAll(regionMeshes, regioni)
    {
        regionMeshes[regioni].setInstance(newInstance);
        regionMeshes[regioni].setPolyFacesBfInstance(newInstance);
    }

    // Write resulting mesh
    Info<< "Writing mesh to " << runTime.name() << nl << endl;
    forAll(regionMeshes, regioni)
    {
        regionMeshes[regioni].write();
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
