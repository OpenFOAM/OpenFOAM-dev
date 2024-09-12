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
    engineCompRatio

Description
    Calculate the compression ratio of the engine combustion chamber

    If the combustion chamber is not the entire mesh a \c cellSet or
    \c cellZone name of the cells in the combustion chamber can be provided.

Usage
    \b engineCompRatio [OPTION]

      - \par -cellSet \<name\>
        Specify the cellSet name of the combustion chamber

      - \par -cellZone zoneName
        Specify the cellZone name of the combustion chamber

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volMesh.H"
#include "cellSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addMeshOption.H"
    #include "addRegionOption.H"

    argList::addOption
    (
        "cellSet",
        "name",
        "Calculate the volume of the specified cellSet"
    );

    argList::addOption
    (
        "cellZone",
        "name",
        "Calculate the volume of the specified cellZone"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSpecifiedMesh.H"

    // Cell labels of the cells in the combustion chamber
    // Defaults to all cells
    labelList ccCells(identityMap(mesh.nCells()));

    word cellSetName;
    word cellZoneName;
    if (args.optionReadIfPresent("cellSet", cellSetName))
    {
        // Read the cellSet for the combustion chamber
        const cellSet ccCellSet(mesh, cellSetName);
        ccCells = ccCellSet.toc();
    }
    else if (args.optionReadIfPresent("cellZone", cellZoneName))
    {
        const cellZone& ccCellZone = mesh.cellZones()[cellZoneName];
        ccCells = ccCellZone;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const scalar eps = 1.0e-10;

    // Search for the CA >= the current CA which is divisible by 180
    int ca0 = -360;
    while (runTime.userTimeValue() > ca0)
    {
        ca0 += 180;
    }
    int ca1 = ca0 + 180;

    if (mag(runTime.userTimeValue() - ca0) > eps)
    {
        const scalar t0 = runTime.userTimeToTime(ca0 - runTime.userTimeValue());
        runTime.setDeltaT(t0);
        runTime++;
        Info<< "Moving to CA = " << runTime.userTimeValue() << endl;
        mesh.move();
    }

    const scalar V0 = gSum(scalarField(mesh.V(), ccCells));
    Info << "Volume at " << runTime.userTimeValue() << "CA = " << V0 << endl;

    if (mag(runTime.userTimeValue()-ca1) > eps)
    {
        const scalar t1 = runTime.userTimeToTime(ca1 - runTime.userTimeValue());
        runTime.setDeltaT(t1);
        runTime++;
        Info<< "Moving to CA = " << runTime.userTimeValue() << endl;
        mesh.move();
    }

    const scalar V1 = gSum(scalarField(mesh.V(), ccCells));
    Info << "Volume at " << runTime.userTimeValue() << "CA = " << V1 << endl;

    const scalar Vmin(min(V0, V1));
    const scalar Vmax(max(V0, V1));

    Info<< "\nCompression ratio = " << Vmax/Vmin << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
