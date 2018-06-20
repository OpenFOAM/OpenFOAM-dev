/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    DPMDyMFoam

Description
    Transient solver for the coupled transport of a single kinematic particle
    cloud including the effect of the volume fraction of particles on the
    continuous phase, with optional mesh motion and mesh topology changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"

#ifdef MPPIC
    #include "basicKinematicMPPICCloud.H"
    #define basicKinematicTypeCloud basicKinematicMPPICCloud
#else
    #include "basicKinematicCollidingCloud.H"
    #define basicKinematicTypeCloud basicKinematicCollidingCloud
#endif

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUcf.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Store the particle positions
        kinematicCloud.storeGlobalPositions();

        mesh.update();

        // Calculate absolute flux from the mapped surface velocity
        phic = mesh.Sf() & Ucf;

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhic.H"
        }

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phic, Uc);

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        continuousPhaseTransport.correct();
        muc = rhoc*continuousPhaseTransport.nu();

        Info<< "Evolving " << kinematicCloud.name() << endl;
        kinematicCloud.evolve();

        // Update continuous phase volume fraction field
        alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaPhic = alphacf*phic;

        fvVectorMatrix cloudSU(kinematicCloud.SU(Uc));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector
            (
                "0",
                cloudSU.dimensions()/dimVolume,
                Zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = Zero;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UcEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                continuousPhaseTurbulence->correct();
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
