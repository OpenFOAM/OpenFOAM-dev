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

Application
    mergeMeshes

Description
    Merges meshes without stitching.

Usage
    \b mergeMeshes [OPTION]

    Options:
      - \par -doc
        Display the documentation in browser

      - \par -srcDoc
        Display the source documentation in browser

      - \par -help
        Print the usage

      - \par -case \<dir\>
        Select a case directory instead of the current working directory

      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -addRegions "'(region1 region2 ... regionN)'"
        Specify list of region meshes to merge.

      - \par -addCases "'(\"casePath1\" \"casePath2\" ... \"casePathN\")'"
        Specify list of case meshes to merge.

      - \par -addCaseRegions "'((\"casePath1\" region1) (\"casePath2\" region2)"
        Specify list of case region meshes to merge.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "mergePolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Merge meshes without stitching");

    argList::noParallel();
    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::addOption
    (
        "addRegions",
        "'(region1 region2 ... regionN)'"
        "list of regions to merge"
    );

    argList::addOption
    (
        "addCases",
        "'(\"casePath1\" \"casePath2\" ... \"casePathN\")'",
        "list of cases to merge"
    );

    argList::addOption
    (
        "addCaseRegions",
        "'((\"casePath1\" region1) (\"casePath2\" region2)"
        "... (\"casePathN\" regionN))'",
        "list of case regions to merge"
    );

    #include "setRootCase.H"

    const wordList regions
    (
        args.optionLookupOrDefault<wordList>("addRegions", wordList::null())
    );

    const fileNameList cases
    (
        args.optionLookupOrDefault<fileNameList>
        (
            "addCases",
            fileNameList::null()
        )
    );

    List<Tuple2<fileName, word>> caseRegions
    (
        args.optionLookupOrDefault<List<Tuple2<fileName, word>>>
        (
            "addCaseRegions",
            List<Tuple2<fileName, word>>::null()
        )
    );

    forAll(cases, i)
    {
        caseRegions.append({cases[i], polyMesh::defaultRegion});
    }

    #include "createTimeNoFunctionObjects.H"

    // Select time if specified
    timeSelector::selectIfPresent(runTime, args);

    #include "createNamedPolyMesh.H"

    const bool overwrite = args.optionFound("overwrite");
    const word oldInstance = mesh.pointsInstance();

    if (!overwrite)
    {
        runTime++;
    }

    // Construct the mergePolyMesh class for the current mesh
    mergePolyMesh mergeMeshes(mesh);


    // Add all the specified region meshes
    forAll(regions, i)
    {
        Info<< "Create polyMesh for region " << regions[i] << endl;
        polyMesh meshToAdd
        (
            IOobject
            (
                regions[i],
                runTime.name(),
                runTime
            )
        );

        Info<< "Adding mesh " << meshToAdd.objectPath() << endl;
        mergeMeshes.addMesh(meshToAdd);
    }

    // Add all the specified case meshes
    forAll(caseRegions, i)
    {
        const fileName& addCase = caseRegions[i].first();
        const word& addRegion = caseRegions[i].second();

        const fileName addCasePath(addCase.path());
        const fileName addCaseName(addCase.name());

        // Construct the time for the new case without reading the controlDict
        Time runTimeToAdd
        (
            Time::controlDictName,
            addCasePath,
            addCaseName,
            false
        );

        Info<< "Create polyMesh for case " << runTimeToAdd.path() << endl;
        polyMesh meshToAdd
        (
            IOobject
            (
                addRegion,
                runTimeToAdd.name(),
                runTimeToAdd
            )
        );

        Info<< "Adding mesh " << meshToAdd.objectPath() << endl;
        mergeMeshes.addMesh(meshToAdd);
    }

    Info << nl << "Merging all meshes" << endl;
    mergeMeshes.merge();

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to " << mesh.facesInstance() << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
