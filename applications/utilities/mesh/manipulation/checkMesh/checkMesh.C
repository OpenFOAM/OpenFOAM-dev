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
    checkMesh

Description
    Checks validity of a mesh.

Usage
    \b checkMesh [OPTION]

    Options:
      - \par noTopology
        Skip checking the mesh topology

      - \par -allTopology
        Check all (including non finite-volume specific) addressing

      - \par -allGeometry
        Check all (including non finite-volume specific) geometry

      - \par -meshQuality
        Check against user defined (in \a system/meshQualityDict) quality
        settings

      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -writeSurfaces
        Reconstruct cellSets and faceSets of problem faces and write to
        postProcessing directory.

      - \par -surfaceFormat <format>
        Format used to write the cellSets and faceSets surfaces
        e.g. vtk or ensight.

      - \par -writeSets
        Reconstruct pointSets of problem points nd write to
        postProcessing directory.

      - \par -setFormat <format>
        Format used to write the pointSets
        e.g. vtk or ensight.

      - \par -nonOrthThreshold <threshold value in degrees>
        Threshold in degrees for reporting non-orthogonality errors,
        default: 70"

      - \par -skewThreshold <threshold value>
        Threshold for reporting skewness errors, default: 4.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "vtkSurfaceWriter.H"
#include "vtkSetWriter.H"
#include "IOdictionary.H"

#include "meshCheck.H"
#include "checkMeshQuality.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "noTopology",
        "skip checking the mesh topology"
    );
    argList::addBoolOption
    (
        "allGeometry",
        "include bounding box checks"
    );
    argList::addBoolOption
    (
        "allTopology",
        "include extra topology checks"
    );
    argList::addBoolOption
    (
        "meshQuality",
        "read user-defined mesh quality criterions from system/meshQualityDict"
    );
    argList::addBoolOption
    (
        "writeSurfaces",
        "reconstruct and write faceSets and cellSets of the problem faces"
    );
    argList::addOption
    (
        "surfaceFormat",
        "surfaceFormat",
        "Format for faceSets and cellSets of the problem faces, defaults to vtk"
    );
    argList::addBoolOption
    (
        "writeSets",
        "reconstruct and write pointSets of the problem points"
    );
    argList::addOption
    (
        "setFormat",
        "setFormat",
        "Format for pointSets of the problem points, defaults to vtk"
    );
    argList::addOption
    (
        "nonOrthThreshold",
        "nonOrthThreshold",
        "Threshold in degrees "
        "for reporting non-orthogonality errors, default: 70"
    );
    argList::addOption
    (
        "skewThreshold",
        "skewThreshold",
        "Threshold for reporting non-orthogonality errors, default: 4"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    const instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createSpecifiedPolyMesh.H"

    const bool noTopology  = args.optionFound("noTopology");
    const bool allGeometry = args.optionFound("allGeometry");
    const bool allTopology = args.optionFound("allTopology");
    const bool meshQuality = args.optionFound("meshQuality");
    const bool writeSurfaces = args.optionFound("writeSurfaces");
    const bool writeSets = args.optionFound("writeSets");

    word surfaceFormat(vtkSurfaceWriter::typeName);
    args.optionReadIfPresent("surfaceFormat", surfaceFormat);

    word setFormat(vtkSetWriter::typeName);
    args.optionReadIfPresent("surfaceFormat", setFormat);

    scalar nonOrthThreshold = 70;
    args.optionReadIfPresent("nonOrthThreshold", nonOrthThreshold);
    nonOrthThreshold = degToRad(nonOrthThreshold);

    scalar skewThreshold = 4;
    args.optionReadIfPresent("skewThreshold", skewThreshold);

    if (noTopology)
    {
        Info<< "Disabling all topology checks." << nl << endl;
    }
    if (allTopology)
    {
        Info<< "Enabling all (cell, face, edge, point) topology checks."
            << nl << endl;
    }
    if (allGeometry)
    {
        Info<< "Enabling all geometry checks." << nl << endl;
    }
    if (meshQuality)
    {
        Info<< "Enabling user-defined geometry checks." << nl << endl;
    }
    if (writeSurfaces)
    {
        Info<< "Reconstructing and writing surface representation of the "
            << "faceSets and cellSets of problem faces in "
            << surfaceFormat << " format" << nl << endl;
    }
    if (writeSets)
    {
        Info<< "Reconstructing and writing the problem points in "
            << setFormat << " format"
            << nl << endl;
    }


    autoPtr<IOdictionary> qualDict;
    if (meshQuality)
    {
        qualDict.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "meshQualityDict",
                    mesh.time().system(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
           )
        );
    }

    autoPtr<Foam::surfaceWriter> surfaceWriter;
    if (writeSurfaces)
    {
        surfaceWriter = surfaceWriter::New
        (
            surfaceFormat,
            mesh.time().writeFormat(),
            mesh.time().writeCompression()
        );
    }

    autoPtr<Foam::setWriter> setWriter;
    if (writeSets)
    {
        setWriter = Foam::setWriter::New
        (
            setFormat,
            mesh.time().writeFormat(),
            mesh.time().writeCompression()
        );
    }


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        polyMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            !timeI
         || state == polyMesh::TOPO_CHANGE
         || state == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "Time = " << runTime.userTimeName() << nl << endl;

            // Reconstruct globalMeshData
            mesh.globalData();

            meshCheck::printMeshStats(mesh, allTopology);

            label nFailedChecks = 0;

            if (!noTopology)
            {
                nFailedChecks += meshCheck::checkTopology
                (
                    mesh,
                    allTopology,
                    surfaceWriter,
                    setWriter
                );
            }

            nFailedChecks += meshCheck::checkGeometry
            (
                mesh,
                allGeometry,
                nonOrthThreshold,
                skewThreshold,
                surfaceWriter,
                setWriter
            );

            if (meshQuality)
            {
                nFailedChecks +=
                    checkMeshQuality(mesh, qualDict(), surfaceWriter);
            }


            // Note: no reduction in nFailedChecks necessary since is
            //       counter of checks, not counter of failed cells,faces etc.

            if (nFailedChecks == 0)
            {
                Info<< "\nMesh OK.\n" << endl;
            }
            else
            {
                Info<< "\nFailed " << nFailedChecks << " mesh checks.\n"
                    << endl;
            }
        }
        else if (state == polyMesh::POINTS_MOVED)
        {
            Info<< "Time = " << runTime.userTimeName() << nl << endl;

            label nFailedChecks = meshCheck::checkGeometry
            (
                mesh,
                allGeometry,
                nonOrthThreshold,
                skewThreshold,
                surfaceWriter,
                setWriter
            );

            if (meshQuality)
            {
                nFailedChecks +=
                    checkMeshQuality(mesh, qualDict(), surfaceWriter);
            }


            if (nFailedChecks)
            {
                Info<< "\nFailed " << nFailedChecks << " mesh checks.\n"
                    << endl;
            }
            else
            {
                Info<< "\nMesh OK.\n" << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
