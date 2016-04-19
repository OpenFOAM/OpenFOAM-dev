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
    pendulum

Description
    Simple swinging pendulum simulation with 1-DoF.  The motion is integrated
    using a symplectic method for just over 2-periods.

\*---------------------------------------------------------------------------*/

#include "rigidBodyModel.H"
#include "masslessBody.H"
#include "rigidBodyModelState.H"
#include "sphere.H"
#include "joints.H"
#include "IFstream.H"

using namespace Foam;
using namespace RBD;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    /*
    bool testMerge = true;

    // Create a model for the pendulum
    rigidBodyModel pendulum;

    // Join a weight to the origin with a centre of mass -1m below the origin
    // by a hinge which rotates about the z-axis
    if (testMerge)
    {
        pendulum.join
        (
            pendulum.bodyID("root"),
            Xt(Zero),
            joint::New(new joints::Rz()),
            autoPtr<rigidBody>(new masslessBody("hinge"))
        );

        pendulum.merge
        (
            pendulum.bodyID("hinge"),
            Xt(vector(0, -1, 0)),
            autoPtr<rigidBody>(new sphere("weight", 1, 0.05))
        );
    }
    else
    {
        pendulum.join
        (
            pendulum.bodyID("root"),
            Xt(Zero),
            joint::New(new joints::Rz()),
            rigidBody::New("pendulum", 1, vector(0, -1, 0), 1e-3*I)
        );
    }
    */

    // Create the pendulum model from dictionary
    rigidBodyModel pendulum(dictionary(IFstream("pendulum")()));

    Info<< pendulum << endl;

    // Create the joint-space state fields
    rigidBodyModelState pendulumState(pendulum);
    scalarField& q = pendulumState.q();
    scalarField& qDot = pendulumState.qDot();
    scalarField& qDdot = pendulumState.qDdot();

    scalarField tau(pendulum.nDoF(), Zero);

    // Set the angle of the pendulum to 0.3rad
    q[0] = 0.3;

    // Set the gravitational acceleration
    pendulum.g() = vector(0, -9.81, 0);

    // Integrate the motion of the pendulum for 4.1s (~2-periods) using a
    // symplectic method
    scalar deltaT = 0.01;
    for (scalar t=0; t<4.1; t+=deltaT)
    {
        qDot += 0.5*deltaT*qDdot;
        q += deltaT*qDot;

        pendulum.forwardDynamics(pendulumState, tau, Field<spatialVector>());

        qDot += 0.5*deltaT*qDdot;

        Info<< "Time " << t << "s, angle = " << q[0] << "rad" << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
