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

      - \par -addMeshes "'(mesh1 mesh2 ... meshN)'"
        Specify list of meshes to merge.

      - \par -addRegions "'(region1 region2 ... regionN)'"
        Specify list of region meshes to merge.

      - \par -addCases "'(\"casePath1\" \"casePath2\" ... \"casePathN\")'"
        Specify list of case meshes to merge.

      - \par -addCaseMeshes "'((\"casePath1\" mesh1) (\"casePath2\" mesh2) ... \
        (\"casePathN\" meshN))'"
        Specify list of case meshes to merge.

      - \par -addCaseRegions "'((\"casePath1\" region1) (\"casePath2\" region2)"
        Specify list of case region meshes to merge.

      - \par -addCaseMeshRegions  "'((\"casePath1\" mesh1 region1) \
        (\"casePath2\" mesh2 region2) ... (\"casePathN\" meshN regionN))'"
        Specify list of case mesh regions to merge.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Tuple3.H"
#include "timeSelector.H"
#include "mergePolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Merge meshes without stitching");

    argList::noParallel();

    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "addOverwriteOption.H"

    timeSelector::addOptions();

    argList::addOption
    (
        "addMeshes",
        "'(mesh1 mesh2 ... meshN)'",
        "list of meshes to merge"
    );

    argList::addOption
    (
        "addRegions",
        "'(region1 region2 ... regionN)'",
        "list of regions to merge"
    );

    argList::addOption
    (
        "addMeshRegions",
        "'((mesh1 region1) (mesh2 region2) ... (mesh3 regionN))'",
        "list of mesh regions to merge"
    );

    argList::addOption
    (
        "addCases",
        "'(\"casePath1\" \"casePath2\" ... \"casePathN\")'",
        "list of cases to merge"
    );

    argList::addOption
    (
        "addCaseMeshes",
        "'((\"casePath1\" mesh1) (\"casePath2\" mesh2)"
        "... (\"casePathN\" meshN))'",
        "list of case meshes to merge"
    );

    argList::addOption
    (
        "addCaseRegions",
        "'((\"casePath1\" region1) (\"casePath2\" region2)"
        "... (\"casePathN\" regionN))'",
        "list of case regions to merge"
    );

    argList::addOption
    (
        "addCaseMeshRegions",
        "'((\"casePath1\" mesh1 region1) (\"casePath2\" mesh2 region2)"
        "... (\"casePathN\" meshN regionN))'",
        "list of case mesh regions to merge"
    );

    #include "setRootCase.H"

    const wordList meshes
    (
        args.optionLookupOrDefault<wordList>
        (
            "addMeshes",
            wordList::null()
        )
    );

    const wordList regions
    (
        args.optionLookupOrDefault<wordList>
        (
            "addRegions",
            wordList::null()
        )
    );

    List<Tuple2<word, word>> meshRegions
    (
        args.optionLookupOrDefault<List<Tuple2<word, word>>>
        (
            "addMeshRegions",
            List<Tuple2<word, word>>::null()
        )
    );

    const fileNameList cases
    (
        args.optionLookupOrDefault<fileNameList>
        (
            "addCases",
            fileNameList::null()
        )
    );

    const List<Tuple2<fileName, word>> caseMeshes
    (
        args.optionLookupOrDefault<List<Tuple2<fileName, word>>>
        (
            "addCaseMeshes",
            List<Tuple2<fileName, word>>::null()
        )
    );

    const List<Tuple2<fileName, word>> caseRegions
    (
        args.optionLookupOrDefault<List<Tuple2<fileName, word>>>
        (
            "addCaseRegions",
            List<Tuple2<fileName, word>>::null()
        )
    );

    List<Tuple3<fileName, word, word>> caseMeshRegions
    (
        args.optionLookupOrDefault<List<Tuple3<fileName, word, word>>>
        (
            "addCaseMeshRegions",
            List<Tuple3<fileName, word, word>>::null()
        )
    );

    forAll(meshes, i)
    {
        meshRegions.append({meshes[i], polyMesh::defaultRegion});
    }

    forAll(regions, i)
    {
        meshRegions.append({word::null, regions[i]});
    }

    forAll(cases, i)
    {
        caseMeshRegions.append({cases[i], word::null, polyMesh::defaultRegion});
    }

    forAll(caseMeshes, i)
    {
        caseMeshRegions.append
        (
            {
                caseMeshes[i].first(),
                caseMeshes[i].second(),
                polyMesh::defaultRegion
            }
        );
    }

    forAll(caseRegions, i)
    {
        caseMeshRegions.append
        (
            {
                caseRegions[i].first(),
                word::null,
                caseRegions[i].second()
            }
        );
    }

    #include "createTimeNoFunctionObjects.H"

    // Select time if specified
    timeSelector::selectIfPresent(runTime, args);

    #include "createSpecifiedPolyMesh.H"

    const bool overwrite = args.optionFound("overwrite");
    const word oldInstance = mesh.pointsInstance();

    if (!overwrite)
    {
        runTime++;
    }

    // Construct the mergePolyMesh class for the current mesh
    mergePolyMesh mergeMeshes(mesh);

    // Add all the specified mesh regions
    forAll(meshRegions, i)
    {
        const fileName addLocal =
            meshRegions[i].first() != word::null
          ? "meshes"/meshRegions[i].first()
          : fileName::null;
        const word& addRegion = meshRegions[i].second();

        Info<< "Reading polyMesh " << addLocal/addRegion << endl;
        polyMesh meshToAdd
        (
            IOobject
            (
                addRegion,
                runTime.name(),
                addLocal,
                runTime
            )
        );

        Info<< "Adding mesh " << meshToAdd.objectPath() << endl;
        mergeMeshes.addMesh(meshToAdd);
    }

    // Add all the specified case meshes
    forAll(caseMeshRegions, i)
    {
        const fileName& addCase = caseMeshRegions[i].first();
        const fileName addLocal =
        caseMeshRegions[i].second() != word::null
          ? "meshes"/caseMeshRegions[i].second()
          : fileName::null;
        const word& addRegion = caseMeshRegions[i].third();

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

        Info<< "Reading polyMesh for case "
            << runTimeToAdd.path()/"constant"/addLocal/addRegion << endl;
        polyMesh meshToAdd
        (
            IOobject
            (
                addRegion,
                runTimeToAdd.name(),
                addLocal,
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
