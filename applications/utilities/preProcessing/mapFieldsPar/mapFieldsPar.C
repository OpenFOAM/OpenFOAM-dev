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
    mapFieldsPar

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMeshToFvMesh.H"
#include "mapGeometricFields.H"
#include "mapClouds.H"
#include "intersectionCellsToCells.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapConsistentMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const word& mapMethod,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Consistently creating and mapping fields for time "
        << srcMesh.time().name() << nl << endl;

    fvMeshToFvMesh interp(srcMesh, tgtMesh, mapMethod);

    Info<< nl << "Mapping geometric fields" << endl;

    mapGeometricFields(interp, wordReList(), selectedFields, noLagrangian);

    if (!noLagrangian)
    {
        mapClouds(interp);
    }
}


void mapMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const HashTable<word>& patchMap,
    const wordReList& cuttingPatches,
    const word& mapMethod,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Creating and mapping fields for time "
        << srcMesh.time().name() << nl << endl;

    fvMeshToFvMesh interp(srcMesh, tgtMesh, mapMethod, patchMap);

    Info<< nl << "Mapping geometric fields" << endl;

    mapGeometricFields(interp, cuttingPatches, selectedFields, noLagrangian);

    if (!noLagrangian)
    {
        mapClouds(interp);
    }
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
        "specify the mapping method"
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

    #include "setRootCase.H"

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

    const word mapMethod
    (
        args.optionLookupOrDefault<word>
        (
            "mapMethod",
            cellsToCellss::intersection::typeName
        )
    );
    Info<< "Mapping method: " << mapMethod << endl;

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    const bool noLagrangian = args.optionFound("noLagrangian");

    #include "createTimes.H"

    HashTable<word> patchMap;
    wordReList cuttingPatches;

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

    fvMesh srcMesh
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.name(),
            runTimeSource
        ),
        false
    );

    srcMesh.postConstruct(false, false, fvMesh::stitchType::nonGeometric);

    fvMesh tgtMesh
    (
        IOobject
        (
            targetRegion,
            runTimeTarget.name(),
            runTimeTarget
        ),
        false
    );

    tgtMesh.postConstruct(false, false, fvMesh::stitchType::nonGeometric);

    Info<< "Source mesh size: "
        << returnReduce(srcMesh.nCells(), sumOp<label>())
        << ", Target mesh size: "
        << returnReduce(tgtMesh.nCells(), sumOp<label>())
        << endl;

    if (consistent)
    {
        mapConsistentMesh
        (
            srcMesh,
            tgtMesh,
            mapMethod,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        mapMesh
        (
            srcMesh,
            tgtMesh,
            patchMap,
            cuttingPatches,
            mapMethod,
            selectedFields,
            noLagrangian
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
