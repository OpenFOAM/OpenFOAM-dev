/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    thermoFoam

Description
    Evolves the thermodynamics on a frozen flow field

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulenceModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "radiationModel.H"
#include "fvIOoptionList.H"
#include "simpleControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createRadiationModel.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEvolving thermodynamics\n" << endl;

    if (mesh.solutionDict().found("SIMPLE"))
    {
        simpleControl simple(mesh);

        while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            while (simple.correctNonOrthogonal())
            {
                #include "EEqn.H"
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            runTime.write();
        }
    }
    else
    {
        pimpleControl pimple(mesh);

        while (runTime.run())
        {
            runTime++;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            while (pimple.correctNonOrthogonal())
            {
                #include "EEqn.H"
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            runTime.write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
