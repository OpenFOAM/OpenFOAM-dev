/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshToMesh.H"
#include "processorPolyPatch.H"
#include "MapMeshes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const word& mapMethod,
    const word& AMIMapMethod,
    const bool subtract,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Consistently creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp(meshSource, meshTarget, mapMethod, AMIMapMethod);

    if (subtract)
    {
        MapMesh<minusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        MapMesh<plusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
}


void mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const word& mapMethod,
    const word& AMIMapMethod,
    const bool subtract,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp
    (
        meshSource,
        meshTarget,
        mapMethod,
        AMIMapMethod,
        patchMap,
        cuttingPatches
    );

    if (subtract)
    {
        MapMesh<minusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        MapMesh<plusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
}


wordList addProcessorPatches
(
    const fvMesh& meshTarget,
    const wordList& cuttingPatches
)
{
    // Add the processor patches to the cutting list
    HashSet<word> cuttingPatchTable;
    forAll(cuttingPatches, i)
    {
        cuttingPatchTable.insert(cuttingPatches[i]);
    }

    const polyBoundaryMesh& pbm = meshTarget.boundaryMesh();

    forAll(pbm, patchI)
    {
        if (isA<processorPolyPatch>(pbm[patchI]))
        {
            const word& patchName = pbm[patchI].name();
            cuttingPatchTable.insert(patchName);
        }
    }

    return cuttingPatchTable.toc();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );

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
        "consistent",
        "source and target geometry and boundary conditions identical"
    );
    argList::addOption
    (
        "mapMethod",
        "word",
        "specify the mapping method (direct|mapNearest|cellVolumeWeight)"
    );
    argList::addOption
    (
        "patchMapMethod",
        "word",
        "specify the patch mapping method (direct|mapNearest|faceAreaWeight)"
    );
    argList::addBoolOption
    (
        "subtract",
        "subtract mapped source from target"
    );
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be mapped. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    argList::addBoolOption
    (
        "noLagrangian",
        "skip mapping lagrangian positions and fields"
    );

    argList args(argc, argv);

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    const fileName casePath = args[1];
    const fileName rootDirSource = casePath.path();
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

    const bool consistent = args.optionFound("consistent");


    word mapMethod = meshToMesh::interpolationMethodNames_
    [
        meshToMesh::imCellVolumeWeight
    ];

    if  (args.optionReadIfPresent("mapMethod", mapMethod))
    {
        Info<< "Mapping method: " << mapMethod << endl;
    }


    word patchMapMethod;
    if (meshToMesh::interpolationMethodNames_.found(mapMethod))
    {
        // Lookup corresponding AMI method
        meshToMesh::interpolationMethod method =
            meshToMesh::interpolationMethodNames_[mapMethod];

        patchMapMethod = AMIPatchToPatchInterpolation::interpolationMethodToWord
        (
            meshToMesh::interpolationMethodAMI(method)
        );
    }

    // Optionally override
    if (args.optionFound("patchMapMethod"))
    {
        patchMapMethod = args["patchMapMethod"];

        Info<< "Patch mapping method: " << patchMapMethod << endl;
    }


    if (patchMapMethod.empty())
    {
        FatalErrorIn(args.executable())
            << "No valid patchMapMethod for method " << mapMethod
            << ". Please supply one through the 'patchMapMethod' option"
            << exit(FatalError);
    }



    const bool subtract = args.optionFound("subtract");
    if (subtract)
    {
        Info<< "Subtracting mapped source field from target" << endl;
    }

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    const bool noLagrangian = args.optionFound("noLagrangian");

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

    #include "setTimeIndex.H"

    Info<< "\nCreate meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            targetRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;

    if (consistent)
    {
        mapConsistentMesh
        (
            meshSource,
            meshTarget,
            mapMethod,
            patchMapMethod,
            subtract,
            selectedFields,
            noLagrangian
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
            mapMethod,
            patchMapMethod,
            subtract,
            selectedFields,
            noLagrangian
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
