/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    denseParticleFoam

Description
    Transient solver for the coupled transport of particle clouds including the
    effect of the volume fraction of particles on the continuous phase, with
    optional mesh motion and mesh topology changes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "viscosityModel.H"
#include "phaseIncompressibleMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fvModels.H"
#include "fvConstraints.H"

#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "fvCorrectPhi.H"
#include "fvcReconstruct.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "parcelClouds.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUcfIfPresent.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Store the particle positions
        clouds.storeGlobalPositions();

        // Move the mesh
        mesh.move();

        if (mesh.changing())
        {
            if (correctPhi)
            {
                #include "correctPhic.H"
            }

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        continuousPhaseViscosity->correct();
        muc = rhoc*continuousPhaseViscosity->nu();

        clouds.evolve();

        // Update continuous phase volume fraction field
        alphac = max(1.0 - clouds.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaPhic = alphacf*phic;

        // Dispersed phase drag force
        volVectorField Fd
        (
            IOobject
            (
                "Fd",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedVector(dimAcceleration, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        // continuous-dispersed phase drag coefficient
        volScalarField Dc
        (
            IOobject
            (
                "Dc",
                runTime.name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        {
            const fvVectorMatrix cloudSU(clouds.SU(Uc));

            // Dispersed phase drag force
            Fd.primitiveFieldRef() = -cloudSU.source()/mesh.V()/rhoc;
            Fd.correctBoundaryConditions();

            // Continuous phase drag coefficient
            Dc.primitiveFieldRef() = -cloudSU.diag()/mesh.V()/rhoc;
            Dc.correctBoundaryConditions();
        }

        // Face continuous-dispersed phase drag coefficient
        const surfaceScalarField Dcf(fvc::interpolate(Dc));

        // Face dispersed phase drag force
        const surfaceScalarField Fdf(fvc::flux(Fd));

        // Effective flux of the dispersed phase
        const surfaceScalarField phid
        (
            Fdf/(Dcf + dimensionedScalar(Dc.dimensions(), small))
        );

        // Face buoyancy force
        const surfaceScalarField Fgf(g & mesh.Sf());

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fvModels.correct();

            #include "UcEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.correctTransport())
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
