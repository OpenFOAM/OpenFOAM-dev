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
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("write", "write mesh/results files");
    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"
    #include "createNamedMesh.H"

    const bool write = args.optionFound("write");
    const bool overwrite = args.optionFound("overwrite");

    if (write || overwrite)
    {
        const word oldInstance = mesh.pointsInstance();

        mesh.setInstance(runTime.name());

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Write resulting mesh
        Info<< "Writing mesh to " << runTime.name() << nl << endl;
        mesh.write();
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
