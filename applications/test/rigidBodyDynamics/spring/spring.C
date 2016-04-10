/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    spring

Description
    Simple weight and damped-spirng simulation with 1-DoF.

\*---------------------------------------------------------------------------*/

#include "rigidBodyModel.H"
#include "masslessBody.H"
#include "sphere.H"
#include "joints.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;
using namespace RBD;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Create the spring model from dictionary
    rigidBodyModel spring(dictionary(IFstream("spring")()));

    Info<< spring << endl;

    // Create the joint-space state fields
    scalarField q(spring.nDoF(), Zero);
    scalarField w(spring.nw(), Zero);
    scalarField qDot(spring.nDoF(), Zero);
    scalarField qDdot(spring.nDoF(), Zero);
    scalarField tau(spring.nDoF(), Zero);
    Field<spatialVector> fx(spring.nBodies(), Zero);

    OFstream qFile("qVsTime");
    OFstream qDotFile("qDotVsTime");

    // Integrate the motion of the spring for 4s using a symplectic method
    scalar deltaT = 0.002;
    for (scalar t=0; t<4; t+=deltaT)
    {
        qDot += 0.5*deltaT*qDdot;
        q += deltaT*qDot;

        // Update the body-state prior to the evaluation of the restraints
        spring.forwardDynamicsCorrection
        (
            q,
            w,
            qDot,
            qDdot
        );

        // Accumulate the restraint forces
        fx = Zero;
        spring.applyRestraints(fx);

        // Calculate the body acceleration for the given state
        // and restraint forces
        spring.forwardDynamics
        (
            q,
            w,
            qDot,
            tau,
            fx,
            qDdot
        );

        // Update the velocity
        qDot += 0.5*deltaT*qDdot;

        // Write the results for graph generation
        // using 'gnuplot spring.gnuplot'
        qFile << t << " " << q[0] << endl;
        qDotFile << t << " " << qDot[0] << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
