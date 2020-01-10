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
    kivaToFoam

Description
    Converts a KIVA3v grid to OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "mergedCyclicPolyPatch.H"
#include "polyMeshUnMergeCyclics.H"
#include "unitConversion.H"

using namespace Foam;

//- Supported KIVA versions
enum kivaVersions
{
    kiva3,
    kiva3v
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption
    (
        "file",
        "name",
        "specify alternative input file name - default is otape17"
    );
    argList::addOption
    (
        "version",
        "version",
        "specify kiva version [kiva3|kiva3v] - default is '3v'"
    );
    argList::addOption
    (
        "zHeadMin",
        "scalar",
        "minimum z-height for transferring liner faces to cylinder-head"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const fileName kivaFileName =
        args.optionLookupOrDefault<fileName>("file", "otape17");

    kivaVersions kivaVersion = kiva3v;
    if (args.optionFound("version"))
    {
        const word versionName = args["version"];

        if (versionName == "kiva3")
        {
            kivaVersion = kiva3;
        }
        else if (versionName == "kiva3v")
        {
            kivaVersion = kiva3v;
        }
        else
        {
            FatalErrorInFunction
                << "KIVA file version " << versionName << " not supported"
                << exit(FatalError);

            args.printUsage();
            FatalError.exit(1);
        }
    }

    scalar zHeadMin = -great;
    args.optionReadIfPresent("zHeadMin", zHeadMin);

    #include "readKivaGrid.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
