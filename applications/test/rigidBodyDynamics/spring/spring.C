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
    Simple weight and damped-spring simulation with 1-DoF.

\*---------------------------------------------------------------------------*/

#include "rigidBodyMotion.H"
#include "masslessBody.H"
#include "sphere.H"
#include "joints.H"
#include "rigidBodyRestraint.H"
#include "rigidBodyModelState.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;
using namespace RBD;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    dictionary springDict(IFstream("spring")());

    // Create the spring model from dictionary
    rigidBodyMotion spring(springDict);

    label nIter(readLabel(springDict.lookup("nIter")));

    Info<< spring << endl;

    // Create the joint-space force field
    scalarField tau(spring.nDoF(), Zero);

    // Create the external body force field
    Field<spatialVector> fx(spring.nBodies(), Zero);

    OFstream qFile("qVsTime");
    OFstream qDotFile("qDotVsTime");

    // Integrate the motion of the spring for 4s
    scalar deltaT = 0.002;
    for (scalar t=0; t<4; t+=deltaT)
    {
        spring.newTime();

        for (label i=0; i<nIter; i++)
        {
            spring.solve(deltaT, tau, fx);
        }

        // Write the results for graph generation
        // using 'gnuplot spring.gnuplot'
        qFile << t << " " << spring.state().q()[0] << endl;
        qDotFile << t << " " << spring.state().qDot()[0] << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
