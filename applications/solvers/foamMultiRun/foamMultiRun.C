/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    foamMultiRun

Description
    Loads and executes an OpenFOAM solver modules for each region of a
    multiregion simulation e.g. for conjugate heat transfer.

    The region solvers are specified in the \c regionSolvers dictionary entry in
    \c controlDict, containing a list of pairs of region and solver names,
    e.g. for a two region case with one fluid region named
    liquid and one solid region named tubeWall:
    \verbatim
        regionSolvers
        {
            liquid          fluid;
            tubeWall        solid;
        }
    \endverbatim

    The \c regionSolvers entry is a dictionary to support name substitutions to
    simplify the specification of a single solver type for a set of
    regions, e.g.
    \verbatim
        fluidSolver     fluid;
        solidSolver     solid;

        regionSolvers
        {
            tube1             $fluidSolver;
            tubeWall1         solid;
            tube2             $fluidSolver;
            tubeWall2         solid;
            tube3             $fluidSolver;
            tubeWall3         solid;
        }
    \endverbatim

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

Usage
    \b foamMultiRun [OPTION]

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To update and run a \c chtMultiRegion case add the following entry to
        the controlDict:
        \verbatim
            regionSolvers
            {
                fluid           fluid;
                solid           solid;
            }
        \endverbatim
        then execute \c foamMultiRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "regionSolvers.H"
#include "pimpleMultiRegionControl.H"
#include "setDeltaT.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Create the region meshes and solvers
    regionSolvers solvers(runTime);

    // Create the outer PIMPLE loop and control structure
    pimpleMultiRegionControl pimple(runTime, solvers);

    // Set the initial time-step
    setDeltaT(runTime, solvers);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        forAll(solvers, i)
        {
            solvers[i].preSolve();
        }

        solvers.setGlobalPrefix();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solvers);

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Multi-region PIMPLE corrector loop
        while (pimple.loop())
        {
            forAll(solvers, i)
            {
                if (solvers[i].pimple.flow())
                {
                    solvers[i].moveMesh();
                }
            }

            forAll(solvers, i)
            {
                if (solvers[i].pimple.flow())
                {
                    solvers[i].motionCorrector();
                }
            }

            forAll(solvers, i)
            {
                if (solvers[i].pimple.models())
                {
                    solvers[i].fvModels().correct();
                }
            }

            forAll(solvers, i)
            {
                solvers[i].prePredictor();
            }

            forAll(solvers, i)
            {
                if
                (
                    solvers[i].pimple.predictTransport()
                 && solvers[i].pimple.flow()
                )
                {
                    solvers[i].momentumTransportPredictor();
                }
            }

            forAll(solvers, i)
            {
                if
                (
                    solvers[i].pimple.predictTransport()
                 && solvers[i].pimple.thermophysics()
                )
                {
                    solvers[i].thermophysicalTransportPredictor();
                }
            }

            forAll(solvers, i)
            {
                if (solvers[i].pimple.flow())
                {
                    solvers[i].momentumPredictor();
                }
            }

            while (pimple.correctEnergy())
            {
                forAll(solvers, i)
                {
                    if (solvers[i].pimple.thermophysics())
                    {
                        solvers[i].thermophysicalPredictor();
                    }
                }
            }

            forAll(solvers, i)
            {
                if (solvers[i].pimple.flow())
                {
                    solvers[i].pressureCorrector();
                }
            }

            forAll(solvers, i)
            {
                if
                (
                    solvers[i].pimple.correctTransport()
                 && solvers[i].pimple.flow()
                )
                {
                    solvers[i].momentumTransportCorrector();
                }
            }

            forAll(solvers, i)
            {
                if
                (
                    solvers[i].pimple.correctTransport()
                 && solvers[i].pimple.thermophysics()
                )
                {
                    solvers[i].thermophysicalTransportCorrector();
                }
            }
        }

        forAll(solvers, i)
        {
            solvers[i].postSolve();
        }

        solvers.setGlobalPrefix();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
