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
    surfaceMeshConvertTesting

Description
    Converts from one surface mesh format to another, but primarily
    used for testing functionality.

Usage
    \b surfaceMeshConvertTesting inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface

      - \par -orient
        Check face orientation on the input surface

      - \par -scale \<scale\>
        Specify a scaling factor for writing the files

      - \par -triSurface
        Use triSurface library for input/output

      - \par -keyed
        Use keyedSurface for input/output

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "surfMesh.H"
#include "surfFields.H"
#include "surfPointFields.H"
#include "PackedBoolList.H"

#include "MeshedSurfaces.H"
#include "UnsortedMeshedSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "convert between surface formats, "
        "but primarily for testing functionality\n"
        "Normally use surfaceMeshConvert instead."
    );

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");

    argList::addBoolOption("clean");
    argList::addBoolOption("orient");
    argList::addBoolOption("surfMesh");
    argList::addBoolOption("triSurface");
    argList::addBoolOption("unsorted");
    argList::addBoolOption("triFace");

    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );

    #include "setRootCase.H"

    const scalar scaleFactor = args.optionLookupOrDefault("scale", 0.0);

    const fileName importName = args[1];
    const fileName exportName = args[2];

    if (importName == exportName)
    {
        FatalErrorInFunction
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }

    if
    (
        !MeshedSurface<face>::canRead(importName, true)
     || !MeshedSurface<face>::canWriteType(exportName.ext(), true)
    )
    {
        return 1;
    }

    if (args.optionFound("triSurface"))
    {
        triSurface surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        if (args.optionFound("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.optionFound("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< endl;
        }

        // write sorted by region
        surf.write(exportName, true);
    }
    else if (args.optionFound("unsorted"))
    {
        UnsortedMeshedSurface<face> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        if (args.optionFound("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.optionFound("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< endl;
        }
        surf.write(exportName);
    }
#if 1
    else if (args.optionFound("triFace"))
    {
        MeshedSurface<triFace> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        if (args.optionFound("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.optionFound("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< endl;
        }
        surf.write(exportName);
    }
#endif
    else
    {
        MeshedSurface<face> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        if (args.optionFound("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.optionFound("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }


        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< endl;
        }
        surf.write(exportName);

        if (args.optionFound("surfMesh"))
        {
            Foam::Time runTime
            (
                args.rootPath(),
                args.caseName()
            );

            // start with "constant"
            runTime.setTime(instant(0, runTime.constant()), 0);

            Info<< "runTime.instance() = " << runTime.instance() << endl;
            Info<< "runTime.timeName() = " << runTime.timeName() << endl;


            Info<< "write MeshedSurface 'yetAnother' via proxy as surfMesh"
                << endl;
            surf.write
            (
                runTime,
                "yetAnother"
            );

            surfMesh surfIn
            (
                IOobject
                (
                    "default",
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );


            MeshedSurface<face> surfIn2(runTime, "foobar");

            Info<<"surfIn2 = " << surfIn2.size() << endl;

            Info<< "surfIn = " << surfIn.size() << endl;


            Info<< "writing surfMesh as obj = oldSurfIn.obj" << endl;

            using Foam::surfMesh;
            surfIn.write(fileName("oldSurfIn.obj"));


            Info<< "runTime.instance() = " << runTime.instance() << endl;

            surfMesh surfOut
            (
                IOobject
                (
                    "mySurf",
                    runTime.instance(),
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                surf.xfer()
            );

            Info<< "writing surfMesh as well: " << surfOut.objectPath() << endl;
            surfOut.write();

            surfLabelField zoneIds
            (
                IOobject
                (
                    "zoneIds",
                    surfOut.instance(),
                    surfOut,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surfOut,
                dimless
            );

            Info<<" surf name= " << surfOut.name() <<nl;
            Info<< "rename to anotherSurf" << endl;
            surfOut.rename("anotherSurf");

            Info<<" surf name= " << surfOut.name() <<nl;

            // advance time to 1
            runTime.setTime(instant(1), 1);
            surfOut.setInstance(runTime.timeName());



            Info<< "writing surfMesh again well: " << surfOut.objectPath()
                << endl;
            surfOut.write();

            // write directly
            surfOut.surfMesh::write(fileName("someName.ofs"));

#if 1
            const surfZoneList& zones = surfOut.surfZones();
            forAll(zones, zoneI)
            {
                SubList<label>
                (
                    zoneIds,
                    zones[zoneI].size(),
                    zones[zoneI].start()
                ) = zoneI;
            }

            Info<< "write zoneIds (for testing only): "
                << zoneIds.objectPath() << endl;
            zoneIds.write();

            surfPointLabelField pointIds
            (
                IOobject
                (
                    "zoneIds.",
//                    "pointIds",
                    surfOut.instance(),
//                    "pointFields",
                    surfOut,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surfOut,
                dimless
            );

            forAll(pointIds, i)
            {
                pointIds[i] = i;
            }

            Info<< "write pointIds (for testing only): "
                << pointIds.objectPath() << endl;
            pointIds.write();

            Info<<"surfMesh with these names: " << surfOut.names() << endl;

#endif
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
