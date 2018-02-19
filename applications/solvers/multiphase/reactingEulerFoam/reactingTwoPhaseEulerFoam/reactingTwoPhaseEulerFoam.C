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
    reactingTwoPhaseEulerFoam

Description
    Solver for a system of 2 compressible fluid phases with a common pressure,
    but otherwise separate properties. The type of phase model is run time
    selectable and can optionally represent multiple species and in-phase
    reactions. The phase system is also run time selectable and can optionally
    represent different types of momentun, heat and mass transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    #include "pUf/createRDeltaTf.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        int nEnergyCorrectors
        (
            pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
        );

        if (LTS)
        {
            #include "setRDeltaT.H"
            if (faceMomentum)
            {
                #include "setRDeltaTf.H"
            }
        }
        else
        {
            #include "CourantNos.H"
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "YEqns.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                #include "EEqns.H"
                #include "pUf/pEqn.H"
            }
            else
            {
                #include "pU/UEqns.H"
                #include "EEqns.H"
                #include "pU/pEqn.H"
            }

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
