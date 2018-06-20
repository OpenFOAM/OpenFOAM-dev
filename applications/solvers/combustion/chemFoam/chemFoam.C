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
    chemFoam

Description
    Solver for chemistry problems, designed for use on single cell cases to
    provide comparison against other chemistry solvers, that uses a single cell
    mesh, and fields created from the initial conditions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiReactionThermo.H"
#include "BasicChemistryModel.H"
#include "reactingMixture.H"
#include "chemistrySolver.H"
#include "OFstream.H"
#include "thermoPhysicsTypes.H"
#include "basicSpecieMixture.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #define CREATE_MESH createSingleCellMesh.H
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createSingleCellMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "readInitialConditions.H"
    #include "createControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"

        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "solveChemistry.H"
        #include "YEqn.H"
        #include "hEqn.H"
        #include "pEqn.H"

        #include "output.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "Number of steps = " << runTime.timeIndex() << endl;
    Info << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
