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
    star4ToFoam

Description
    Converts a Star-CD (v4) pro-STAR mesh into OpenFOAM format.

Usage
    \b star4ToFoam [OPTION] ccmMesh

    Options:
      - \par -ascii
        Write in ASCII format instead of binary

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 0.001 (scale \em [mm] to \em [m]).

      - \par -solids
        Treat any solid cells present just like fluid cells.
        The default is to discard them.

Note
    Baffles are written as interfaces for later use

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "STARCDMeshReader.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "convert pro-STAR (v4) mesh to OpenFOAM"
    );

    argList::noParallel();
    argList::validArgs.append("pro-STAR prefix");
    argList::addBoolOption
    (
        "ascii",
        "write in ASCII instead of binary format"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 0.001 ([mm] to [m])"
    );
    argList::addBoolOption
    (
        "solids",
        "retain solid cells and treat them like fluid cells"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    // default rescale from [mm] to [m]
    scalar scaleFactor = args.optionLookupOrDefault("scale", 0.001);
    if (scaleFactor <= 0)
    {
        scaleFactor = 1;
    }

    meshReaders::STARCD::keepSolids = args.optionFound("solids");

    // default to binary output, unless otherwise specified
    IOstream::streamFormat format = IOstream::BINARY;
    if (args.optionFound("ascii"))
    {
        format = IOstream::ASCII;
    }

    // increase the precision of the points data
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    // remove extensions and/or trailing '.'
    const fileName prefix = fileName(args[1]).lessExt();

    meshReaders::STARCD reader(prefix, runTime, scaleFactor);

    autoPtr<polyMesh> mesh = reader.mesh(runTime);
    reader.writeMesh(mesh, format);


    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
