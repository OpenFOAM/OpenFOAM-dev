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
    collapseEdges

Description
    Collapses short edges and combines edges that are in line.

    - collapse short edges. Length of edges to collapse provided as argument.
    - merge two edges if they are in line. Maximum angle provided as argument.
    - remove unused points.
    - collapse faces:
        - with small areas to a single point
        - that have a high aspect ratio (i.e. sliver face) to a single edge

    Optionally checks the resulting mesh for bad faces and reduces the desired
    face length factor for those faces attached to the bad faces.

    When collapsing an edge with one point on the boundary it will leave
    the boundary point intact. When both points inside it chooses random. When
    both points on boundary random again.

Usage
    - collapseEdges [OPTION]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "polyTopoChange.H"
#include "fvMesh.H"
#include "polyMeshFilter.H"
#include "faceSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, false);
    argList::addNote
    (
        "Collapses small edges to a point.\n"
        "Optionally collapse small faces to a point and thin faces to an edge."
    );

    argList::addBoolOption
    (
        "collapseFaces",
        "Collapse small and sliver faces as well as small edges"
    );

    argList::addOption
    (
        "collapseFaceSet",
        "faceSet",
        "Collapse faces that are in the supplied face set"
    );

    #include "addDictOption.H"
    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const word dictName("collapseDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    IOdictionary collapseDict(dictIO);

    const bool overwrite = args.optionFound("overwrite");

    const bool collapseFaces = args.optionFound("collapseFaces");
    const bool collapseFaceSet = args.optionFound("collapseFaceSet");

    if (collapseFaces && collapseFaceSet)
    {
        FatalErrorInFunction
            << "Both face zone collapsing and face collapsing have been"
            << "selected. Choose only one of:" << nl
            << "    -collapseFaces" << nl
            << "    -collapseFaceSet <faceSet>"
            << abort(FatalError);
    }


    // maintain indirectPatchFaces if it is there (default) or force
    // (if collapseFaceSet option provided)
    word faceSetName("indirectPatchFaces");
    IOobject::readOption readFlag = IOobject::READ_IF_PRESENT;

    if (args.optionReadIfPresent("collapseFaceSet", faceSetName))
    {
        readFlag = IOobject::MUST_READ;
    }



    labelIOList pointPriority
    (
        IOobject
        (
            "pointPriority",
            runTime.timeName(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nPoints(), labelMin)
    );
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        autoPtr<polyMeshFilter> meshFilterPtr;

        label nBadFaces = 0;

        faceSet indirectPatchFaces
        (
            mesh,
            faceSetName,
            readFlag,
            IOobject::AUTO_WRITE
        );
        Info<< "Read faceSet " << indirectPatchFaces.name()
            << " with "
            << returnReduce(indirectPatchFaces.size(), sumOp<label>())
            << " faces" << endl;


        {
            meshFilterPtr.set
            (
                new polyMeshFilter(mesh, pointPriority, collapseDict)
            );
            polyMeshFilter& meshFilter = meshFilterPtr();

            // newMesh will be empty until it is filtered
            const autoPtr<fvMesh>& newMesh = meshFilter.filteredMesh();

            // Filter small edges only. This reduces the number of faces so that
            // the face filtering is sped up.
            nBadFaces = meshFilter.filterEdges(0);
            {
                polyTopoChange meshMod(newMesh());

                meshMod.changeMesh(mesh, false);

                polyMeshFilter::copySets(newMesh(), mesh);
            }

            pointPriority = meshFilter.pointPriority();
        }

        if (collapseFaceSet)
        {
            meshFilterPtr.reset
            (
                new polyMeshFilter(mesh, pointPriority, collapseDict)
            );
            polyMeshFilter& meshFilter = meshFilterPtr();

            const autoPtr<fvMesh>& newMesh = meshFilter.filteredMesh();

            // Filter faces. Pass in the number of bad faces that are present
            // from the previous edge filtering to use as a stopping criterion.
            meshFilter.filter(indirectPatchFaces);
            {
                polyTopoChange meshMod(newMesh);

                meshMod.changeMesh(mesh, false);

                polyMeshFilter::copySets(newMesh(), mesh);
            }

            pointPriority = meshFilter.pointPriority();
        }

        if (collapseFaces)
        {
            meshFilterPtr.reset
            (
                new polyMeshFilter(mesh, pointPriority, collapseDict)
            );
            polyMeshFilter& meshFilter = meshFilterPtr();

            const autoPtr<fvMesh>& newMesh = meshFilter.filteredMesh();

            // Filter faces. Pass in the number of bad faces that are present
            // from the previous edge filtering to use as a stopping criterion.
            meshFilter.filter(nBadFaces);
            {
                polyTopoChange meshMod(newMesh);

                meshMod.changeMesh(mesh, false);

                polyMeshFilter::copySets(newMesh(), mesh);
            }

            pointPriority = meshFilter.pointPriority();
        }

        // Write resulting mesh
        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        Info<< nl << "Writing collapsed mesh to time "
            << runTime.timeName() << nl << endl;

        mesh.write();
        pointPriority.write();
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
