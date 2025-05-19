/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    createZones

Description
    Executes the zoneGenerators specified and writes the zones generated

Usage
    \b createZones [OPTION]

    Options:
      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -mesh \<name\>
        Specify an alternative mesh.

      - \par -dict \<filename\>
        Specify alternative dictionary for the zoneGenerators,
        defaults to system/createZonesDict

      - \par -clear
        Clear all registered zones before generating new zones.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "zoneGeneratorList.H"
#include "zonesGenerator.H"
#include "systemDict.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "clear",
        "clear all existing zones"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createSpecifiedPolyMesh.H"

    if (args.optionFound("clear"))
    {
        mesh.pointZones().clear();
        mesh.cellZones().clear();
        mesh.faceZones().clear();
    }

    const dictionary zoneGeneratorsDict
    (
        systemDict("createZonesDict", args, mesh)
    );

    zoneGeneratorList zoneGenerators(mesh, zoneGeneratorsDict);

    // Generate and register the zones
    zoneGenerators.generate();

    // Generate the zones specified in constant/zonesGenerator
    // zonesGenerator::New(mesh);

    mesh.pointZones().write();
    mesh.cellZones().write();
    mesh.faceZones().write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
