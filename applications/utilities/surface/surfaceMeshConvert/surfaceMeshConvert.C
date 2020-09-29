/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    surfaceMeshConvert

Description
    Converts between surface formats with optional scaling or
    transformations (rotate/translate) on a coordinateSystem.

Usage
    \b surfaceMeshConvert inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface.

      - \par -scaleIn \<scale\>
        Specify a scaling factor when reading files.

      - \par -scaleOut \<scale\>
        Specify a scaling factor when writing files.

      - \par -from \<coordinateSystem\>
        Specify a coordinate System when reading files.

      - \par -to \<coordinateSystem\>
        Specify a coordinate System when writing files.

      - \par -tri
        Triangulate surface.

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "MeshedSurfaces.H"
#include "coordinateSystems.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "convert between surface formats"
    );

    argList::noParallel();
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");

    argList::addBoolOption
    (
        "clean",
        "perform some surface checking/cleanup on the input surface"
    );
    argList::addOption
    (
        "scaleIn",
        "factor",
        "geometry scaling factor on input"
    );
    argList::addOption
    (
        "scaleOut",
        "factor",
        "geometry scaling factor on output"
    );
    argList::addOption
    (
        "from",
        "system",
        "specify the source coordinate system, applied after '-scaleIn'"
    );
    argList::addOption
    (
        "to",
        "system",
        "specify the target coordinate system, applied before '-scaleOut'"
    );
    argList::addBoolOption
    (
        "tri",
        "triangulate surface"
    );


    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const fileName importName = args[1];
    const fileName exportName = args[2];

    // disable inplace editing
    if (importName == exportName)
    {
        FatalErrorInFunction
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }

    // check that reading/writing is supported
    if
    (
        !MeshedSurface<face>::canRead(importName, true)
     || !MeshedSurface<face>::canWriteType(exportName.ext(), true)
    )
    {
        return 1;
    }


    // Get the coordinate transformations
    autoPtr<coordinateSystem> fromCsys;
    autoPtr<coordinateSystem> toCsys;

    if (args.optionFound("from") || args.optionFound("to"))
    {
        coordinateSystems::coordinateSystems csLst(runTime);

        if (args.optionFound("from"))
        {
            const word csName = args["from"];
            fromCsys = csLst[csName].clone();
        }

        if (args.optionFound("to"))
        {
            const word csName = args["to"];
            toCsys = csLst[csName].clone();
        }
    }


    {
        MeshedSurface<face> surf(importName);

        if (args.optionFound("clean"))
        {
            surf.cleanup(true);
        }

        scalar scaleIn = 0;
        if (args.optionReadIfPresent("scaleIn", scaleIn) && scaleIn > 0)
        {
            Info<< " -scaleIn " << scaleIn << endl;
            surf.scalePoints(scaleIn);
        }


        if (fromCsys.valid())
        {
            Info<< " -from " << fromCsys().name() << endl;
            tmp<pointField> tpf = fromCsys().localPosition(surf.points());
            surf.movePoints(tpf());
        }

        if (toCsys.valid())
        {
            Info<< " -to " << toCsys().name() << endl;
            tmp<pointField> tpf = toCsys().globalPosition(surf.points());
            surf.movePoints(tpf());
        }

        scalar scaleOut = 0;
        if (args.optionReadIfPresent("scaleOut", scaleOut) && scaleOut > 0)
        {
            Info<< " -scaleOut " << scaleOut << endl;
            surf.scalePoints(scaleOut);
        }

        if (args.optionFound("tri"))
        {
            Info<< "triangulate" << endl;
            surf.triangulate();
        }

        Info<< "writing " << exportName;
        surf.write(exportName);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
