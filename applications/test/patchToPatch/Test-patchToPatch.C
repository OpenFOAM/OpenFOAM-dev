/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "cpuTime.H"
#include "patchToPatch.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;


int main(int argc, char *argv[])
{
    argList::validArgs.append("source");
    argList::validArgs.append("target");
    argList::validArgs.append("method");
    argList::validArgs.append("reverse");

    argList::addOption
    (
        "sourceCase",
        "dir",
        "The directory of the case with the source patch"
    );
    argList::addOption
    (
        "sourceRegion",
        "name",
        "The region with the source patch"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    // Optionally read a different mesh for the source
    autoPtr<Time> srcRunTimePtr;
    autoPtr<polyMesh> srcMeshPtr;
    if (args.optionFound("sourceCase") || args.optionFound("sourceRegion"))
    {
        const string tgtCase = getEnv("FOAM_CASE");
        const string tgtCaseName = getEnv("FOAM_CASENAME");

        fileName sourceCase =
            args.optionLookupOrDefault<fileName>("sourceCase", tgtCase);
        sourceCase.clean();
        const fileName sourceCaseName =
            Pstream::parRun()
          ? fileName(sourceCase.name())/args.caseName().name()
          : fileName(sourceCase.name());

        setEnv("FOAM_CASE", sourceCase, true);
        setEnv("FOAM_CASENAME", sourceCase.name(), true);

        srcRunTimePtr.set
        (
            new Time(sourceCase.path(), sourceCaseName)
        );

        setEnv("FOAM_CASE", tgtCase, true);
        setEnv("FOAM_CASENAME", tgtCaseName, true);

        srcMeshPtr.set
        (
            new polyMesh
            (
                Foam::IOobject
                (
                    args.optionLookupOrDefault<word>
                    (
                        "sourceRegion",
                        Foam::polyMesh::defaultRegion
                    ),
                    srcRunTimePtr->name(),
                    srcRunTimePtr(),
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }

    const polyMesh& srcMesh = srcMeshPtr.valid() ? srcMeshPtr() : mesh;
    const polyMesh& tgtMesh = mesh;

    const polyPatch& srcPatch = srcMesh.boundaryMesh()[args[1]];
    const polyPatch& tgtPatch = tgtMesh.boundaryMesh()[args[2]];
    const word& method = args[3];
    const bool reverse = args.argRead<bool>(4);

    cpuTime time;

    patchToPatch::New(method, reverse)->update
    (
        srcPatch,
        srcPatch.pointNormals(),
        tgtPatch
    );

    Info<< nl << patchToPatch::typeName << ": Completed in "
        << time.cpuTimeIncrement() << " s" << nl << endl;

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
