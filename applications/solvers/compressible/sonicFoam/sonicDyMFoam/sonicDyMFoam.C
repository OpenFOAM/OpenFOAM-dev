/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    sonicDyMFoam

Group
    grpCompressibleSolvers grpMovingMeshSolvers

Description
    Transient solver for trans-sonic/supersonic, turbulent flow of a
    compressible gas, with optional mesh motion and mesh topology changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUf.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        {
            // Store divrhoU from the previous mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phi has the
            // same divergence
            volScalarField divrhoU
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );

            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"

            runTime++;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Store momentum to set rhoUf for introduced faces.
            volVectorField rhoU("rhoU", rho*U);

            // Do any mesh changes
            mesh.update();

            if (mesh.changing())
            {
                MRF.update();

                if (correctPhi)
                {
                    // Calculate absolute flux from the mapped surface velocity
                    phi = mesh.Sf() & rhoUf;

                    #include "correctPhi.H"

                    // Make the fluxes relative to the mesh-motion
                    fvc::makeRelative(phi, rho, U);
                }
            }

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        #include "rhoEqn.H"
        Info<< "rhoEqn max/min : " << max(rho).value()
            << " " << min(rho).value() << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
