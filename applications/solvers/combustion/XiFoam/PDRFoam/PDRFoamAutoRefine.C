/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    PDRFoam

Description
    Solver for compressible premixed/partially-premixed combustion with
    turbulence modelling.

    Combusting RANS code using the b-Xi two-equation model.
    Xi may be obtained by either the solution of the Xi transport
    equation or from an algebraic expression.  Both approaches are
    based on Gulder's flame speed correlation which has been shown
    to be appropriate by comparison with the results from the
    spectral model.

    Strain effects are incorporated directly into the Xi equation
    but not in the algebraic approximation.  Further work need to be
    done on this issue, particularly regarding the enhanced removal rate
    caused by flame compression.  Analysis using results of the spectral
    model will be required.

    For cases involving very lean Propane flames or other flames which are
    very strain-sensitive, a transport equation for the laminar flame
    speed is present.  This equation is derived using heuristic arguments
    involving the strain time scale and the strain-rate at extinction.
    the transport velocity is the same as that for the Xi equation.

    For large flames e.g. explosions additional modelling for the flame
    wrinkling due to surface instabilities may be applied.

    PDR (porosity/distributed resistance) modelling is included to handle
    regions containing blockages which cannot be resolved by the mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiuMulticomponentThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermoThermophysicalTransportModel.H"
#include "laminarFlameSpeed.H"
#include "XiModel.H"
#include "PDRDragModel.H"
#include "ignition.H"
#include "Switch.H"
#include "bound.H"
#include "pimpleControl.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "readCombustionProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();
    scalar StCoNum = 0.0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    bool hasChanged = false;

    while (pimple.run(runTime))
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        // Indicators for refinement. Note: before runTime++
        // only for postprocessing reasons.
        tmp<volScalarField> tmagGradP = mag(fvc::grad(p));
        volScalarField normalisedGradP
        (
            "normalisedGradP",
            tmagGradP()/max(tmagGradP())
        );
        normalisedGradP.writeOpt() = IOobject::AUTO_WRITE;
        tmagGradP.clear();

        runTime++;

        Info<< "\n\nTime = " << runTime.name() << endl;

        {
            // Make the fluxes absolute
            fvc::makeAbsolute(phi, rho, U);

            // Test : disable refinement for some cells
            // PackedBoolList& protectedCell =
            //     refCast<dynamicRefineFvMesh>(mesh).protectedCell();

            // if (protectedCell.empty())
            // {
            //     protectedCell.setSize(mesh.nCells());
            //     protectedCell = 0;
            // }

            // forAll(betav, celli)
            // {
            //     if (betav[celli] < 0.99)
            //     {
            //         protectedCell[celli] = 1;
            //     }
            // }

            // Flux estimate for introduced faces.
            volVectorField rhoU("rhoU", rho*U);

            // Do any mesh changes
            bool meshChanged = mesh.update();


            if (meshChanged)
            {
                hasChanged = true;
            }

            if (runTime.write() && hasChanged)
            {
                betav.write();
                Lobs.write();
                CT.write();
                drag->writeFields();
                flameWrinkling->writeFields();
                hasChanged = false;
            }

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);
        }


        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fvModels.correct();

            if (pimple.predictTransport())
            {
                turbulence->predict();
                thermophysicalTransport->predict();
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "bEqn.H"
                #include "ftEqn.H"
                #include "EauEqn.H"
                #include "EaEqn.H"

                if (!ign.ignited())
                {
                    thermo.heu() == thermo.he();
                }

                #include "pEqn.H"
            }

            if (pimple.correctTransport())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        runTime.write();

        Info<< "\nExecutionTime = "
             << runTime.elapsedCpuTime()
             << " s\n" << endl;
    }

    Info<< "\n end\n";

    return 0;
}


// ************************************************************************* //
