/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    mapFields

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.
    Parallel and non-parallel cases are handled without the need to reconstruct
    them first.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "mapMeshes.H"
#include "decompositionMethod.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );
    argList::noParallel();
    argList::validArgs.append("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "specify the source time"
    );
    argList::addOption
    (
        "sourceRegion",
        "word",
        "specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "specify the target region"
    );
    argList::addBoolOption
    (
        "parallelSource",
        "the source is decomposed"
    );
    argList::addBoolOption
    (
        "parallelTarget",
        "the target is decomposed"
    );
    argList::addBoolOption
    (
        "consistent",
        "source and target geometry and boundary conditions identical"
    );
    argList::addOption
    (
        "mapMethod",
        "word",
        "specify the mapping method"
    );

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
    word sourceRegion = fvMesh::defaultRegion;
    if (args.optionFound("sourceRegion"))
    {
        sourceRegion = args["sourceRegion"];
        Info<< "Source region: " << sourceRegion << endl;
    }

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = fvMesh::defaultRegion;
    if (args.optionFound("targetRegion"))
    {
        targetRegion = args["targetRegion"];
        Info<< "Target region: " << targetRegion << endl;
    }

    const bool parallelSource = args.optionFound("parallelSource");
    const bool parallelTarget = args.optionFound("parallelTarget");
    const bool consistent = args.optionFound("consistent");

    meshToMesh0::order mapOrder = meshToMesh0::INTERPOLATE;
    if (args.optionFound("mapMethod"))
    {
        const word mapMethod(args["mapMethod"]);
        if (mapMethod == "mapNearest")
        {
            mapOrder = meshToMesh0::MAP;
        }
        else if (mapMethod == "interpolate")
        {
            mapOrder = meshToMesh0::INTERPOLATE;
        }
        else if (mapMethod == "cellPointInterpolate")
        {
            mapOrder = meshToMesh0::CELL_POINT_INTERPOLATE;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mapMethod " << mapMethod << ". Valid options are: "
                << "mapNearest, interpolate and cellPointInterpolate"
                << exit(FatalError);
        }

        Info<< "Mapping method: " << mapMethod << endl;
    }


    #include "createTimes.H"

    HashTable<word> patchMap;
    wordList cuttingPatches;

    if (!consistent)
    {
        IOdictionary mapFieldsDict
        (
            IOobject
            (
                "mapFieldsDict",
                runTimeTarget.system(),
                runTimeTarget,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        mapFieldsDict.lookup("patchMap") >> patchMap;
        mapFieldsDict.lookup("cuttingPatches") >>  cuttingPatches;
    }

    if (parallelSource && !parallelTarget)
    {
        const int nProcs
        (
            decompositionMethod::decomposeParDict(runTimeSource).lookup<int>
            (
                "numberOfSubdomains"
            )
        );

        Info<< "Create target mesh\n" << endl;

        fvMesh meshTarget
        (
            IOobject
            (
                targetRegion,
                runTimeTarget.name(),
                runTimeTarget
            ),
            false
        );

        Info<< "Target mesh size: " << meshTarget.nCells() << endl;

        for (int proci=0; proci<nProcs; proci++)
        {
            Info<< nl << "Source processor " << proci << endl;

            Time runTimeSource
            (
                Time::controlDictName,
                rootDirSource,
                caseDirSource/fileName(word("processor") + name(proci))
            );

            #include "setTimeIndex.H"

            fvMesh meshSource
            (
                IOobject
                (
                    sourceRegion,
                    runTimeSource.name(),
                    runTimeSource
                ),
                false
            );

            Info<< "mesh size: " << meshSource.nCells() << endl;

            if (consistent)
            {
                mapConsistentSubMesh
                (
                    meshSource,
                    meshTarget,
                    mapOrder
                );
            }
            else
            {
                mapSubMesh
                (
                    meshSource,
                    meshTarget,
                    patchMap,
                    cuttingPatches,
                    mapOrder
                );
            }
        }
    }
    else if (!parallelSource && parallelTarget)
    {
        const int nProcs
        (
            decompositionMethod::decomposeParDict(runTimeSource).lookup<int>
            (
                "numberOfSubdomains"
            )
        );

        Info<< "Create source mesh\n" << endl;

        #include "setTimeIndex.H"

        fvMesh meshSource
        (
            IOobject
            (
                sourceRegion,
                runTimeSource.name(),
                runTimeSource
            ),
            false
        );

        Info<< "Source mesh size: " << meshSource.nCells() << endl;

        for (int proci=0; proci<nProcs; proci++)
        {
            Info<< nl << "Target processor " << proci << endl;

            Time runTimeTarget
            (
                Time::controlDictName,
                rootDirTarget,
                caseDirTarget/fileName(word("processor") + name(proci))
            );

            fvMesh meshTarget
            (
                IOobject
                (
                    targetRegion,
                    runTimeTarget.name(),
                    runTimeTarget
                ),
                false
            );

            Info<< "mesh size: " << meshTarget.nCells() << endl;

            if (consistent)
            {
                mapConsistentSubMesh
                (
                    meshSource,
                    meshTarget,
                    mapOrder
                );
            }
            else
            {
                mapSubMesh
                (
                    meshSource,
                    meshTarget,
                    patchMap,
                    addProcessorPatches(meshTarget, cuttingPatches),
                    mapOrder
                );
            }
        }
    }
    else if (parallelSource && parallelTarget)
    {
        const int nProcsSource
        (
            decompositionMethod::decomposeParDict(runTimeSource).lookup<int>
            (
                "numberOfSubdomains"
            )
        );

        const int nProcsTarget
        (
            decompositionMethod::decomposeParDict(runTimeTarget).lookup<int>
            (
                "numberOfSubdomains"
            )
        );

        List<boundBox> bbsTarget(nProcsTarget);
        List<bool> bbsTargetSet(nProcsTarget, false);

        for (int procISource=0; procISource<nProcsSource; procISource++)
        {
            Info<< nl << "Source processor " << procISource << endl;

            Time runTimeSource
            (
                Time::controlDictName,
                rootDirSource,
                caseDirSource/fileName(word("processor") + name(procISource))
            );

            #include "setTimeIndex.H"

            fvMesh meshSource
            (
                IOobject
                (
                    sourceRegion,
                    runTimeSource.name(),
                    runTimeSource
                ),
                false
            );

            Info<< "mesh size: " << meshSource.nCells() << endl;

            boundBox bbSource(meshSource.bounds());

            for (int procITarget=0; procITarget<nProcsTarget; procITarget++)
            {
                if
                (
                    !bbsTargetSet[procITarget]
                  || (
                      bbsTargetSet[procITarget]
                   && bbsTarget[procITarget].overlaps(bbSource)
                     )
                )
                {
                    Info<< nl << "Target processor " << procITarget << endl;

                    Time runTimeTarget
                    (
                        Time::controlDictName,
                        rootDirTarget,
                        caseDirTarget/fileName(word("processor")
                      + name(procITarget))
                    );

                    fvMesh meshTarget
                    (
                        IOobject
                        (
                            targetRegion,
                            runTimeTarget.name(),
                            runTimeTarget
                        ),
                        false
                    );

                    Info<< "mesh size: " << meshTarget.nCells() << endl;

                    bbsTarget[procITarget] = meshTarget.bounds();
                    bbsTargetSet[procITarget] = true;

                    if (bbsTarget[procITarget].overlaps(bbSource))
                    {
                        if (consistent)
                        {
                            mapConsistentSubMesh
                            (
                                meshSource,
                                meshTarget,
                                mapOrder
                            );
                        }
                        else
                        {
                            mapSubMesh
                            (
                                meshSource,
                                meshTarget,
                                patchMap,
                                addProcessorPatches(meshTarget, cuttingPatches),
                                mapOrder
                            );
                        }
                    }
                }
            }
        }
    }
    else
    {
        #include "setTimeIndex.H"

        Info<< "Create meshes\n" << endl;

        fvMesh meshSource
        (
            IOobject
            (
                sourceRegion,
                runTimeSource.name(),
                runTimeSource
            ),
            false
        );

        fvMesh meshTarget
        (
            IOobject
            (
                targetRegion,
                runTimeTarget.name(),
                runTimeTarget
            ),
            false
        );

        Info<< "Source mesh size: " << meshSource.nCells() << tab
            << "Target mesh size: " << meshTarget.nCells() << nl << endl;

        if (consistent)
        {
            mapConsistentMesh(meshSource, meshTarget, mapOrder);
        }
        else
        {
            mapSubMesh
            (
                meshSource,
                meshTarget,
                patchMap,
                cuttingPatches,
                mapOrder
            );
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
