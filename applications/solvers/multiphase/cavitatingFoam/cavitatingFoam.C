/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    cavitatingFoam

Description
    Transient cavitation code based on the homogeneous equilibrium model
    from which the compressibility of the liquid/vapour "mixture" is obtained,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "cavitatingTwoPhaseMixture.H"
#include "incompressibleMomentumTransportModels.H"
#include "fvCorrectPhi.H"
#include "pimpleControl.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "constrainHbyA.H"

#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "createPcorrTypes.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readControls.H"

        {
            #include "CourantNo.H"
            #include "setDeltaT.H"

            // Update the mesh for topology change, mesh to mesh mapping
            mesh.update();

            runTime++;

            Info<< "Time = " << runTime.userTimeName() << nl << endl;

            // Move the mesh
            mesh.move();

            if (mesh.changing() && correctPhi)
            {
                #include "correctPhi.H"
            }
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "rhoEqn.H"

            mixture.correct();

            if (pimple.predictTransport())
            {
                turbulence->predict();
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.correctTransport())
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
