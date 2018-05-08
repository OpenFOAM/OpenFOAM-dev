/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    RBD

Description
    Does an RBD simulation and outputs the result as a gnuplot animation

\*---------------------------------------------------------------------------*/

#include "rigidBodyMotion.H"
#include "IFstream.H"
#include "OFstream.H"
#include "boundBox.H"
#include "argList.H"

using namespace Foam;
using namespace RBD;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Create the input dictionary
    argList::validArgs.append("dictionary");
    argList args(argc, argv);
    const word dictName(args[1]);
    Info<< "Reading " << dictName << nl << endl;
    const dictionary dict = IFstream(dictName)();

    // Read the model, time controls and plot outlines from the dictionary
    rigidBodyMotion motion(dict);
    const scalar deltaT(readScalar(dict.lookup("deltaT")));
    const label nIter(dict.lookupOrDefault<label>("nIter", 1));
    const scalar endTime(readScalar(dict.lookup("endTime")));
    const dictionary& bodiesDict = dict.subDict("bodies");
    List<vectorField> outlines;
    labelList outlineBodyIDs;
    forAll(motion.bodies(), bodyID)
    {
        const word& bodyName = motion.bodies()[bodyID].name();
        if (bodiesDict.isDict(bodyName))
        {
            const dictionary& bodyDict = bodiesDict.subDict(bodyName);
            if (bodyDict.found("outline"))
            {
                outlines.append(vectorField(bodyDict.lookup("outline")));
                outlineBodyIDs.append(bodyID);
            }
        }
    }

    // Set up motion fields
    const label nDoF = motion.nDoF();
    Info << nDoF << " degrees of freedom" << endl;
    scalarField tau(nDoF, Zero);
    Field<spatialVector> fx(motion.nBodies(), Zero);

    // Initialise output files and the bound box
    OFstream dataFile(dictName + ".dat");
    OFstream animationFile(dictName + ".animate");
    boundBox animationBox;

    dataFile<< "# time(s)";
    forAll(motion.state().q(), stateI)
    {
        dataFile << " q_" << stateI;
    }
    forAll(motion.state().q(), stateI)
    {
        dataFile << " qDot_" << stateI;
    }
    forAll(motion.state().q(), stateI)
    {
        dataFile << " qDdot_" << stateI;
    }
    dataFile<< endl;
    animationFile
        << "#!/usr/bin/env gnuplot" << endl << endl
        << "$data << end" << endl;

    // Run the RBD simulation
    for (scalar t = 0; t <= endTime + 0.5*deltaT; t += deltaT)
    {
        Info().stdStream() << "\33[2KTime = " << t << '\r' << std::flush;

        motion.newTime();
        for (label i = 0; i < nIter; ++ i)
        {
            motion.solve(t + deltaT, deltaT, tau, fx);
        }

        // Write the current state
        dataFile<< t;
        forAll(motion.state().q(), stateI)
        {
            dataFile << ' ' << motion.state().q()[stateI];
        }
        forAll(motion.state().q(), stateI)
        {
            dataFile << ' ' << motion.state().qDot()[stateI];
        }
        forAll(motion.state().q(), stateI)
        {
            dataFile << ' ' << motion.state().qDdot()[stateI];
        }
        dataFile<< endl;

        // Write the bodies' outlines at the current transformation
        forAll(outlines, outlineI)
        {
            const label bodyID = outlineBodyIDs[outlineI];
            const spatialTransform& invX0 = motion.X0(bodyID).inv();
            forAll(outlines[outlineI], i)
            {
                const vector p = invX0.transformPoint(outlines[outlineI][i]);
                animationFile<< t << ' ' << p.x() << ' ' << p.y() << endl;
                animationBox.min() = min(animationBox.min(), p);
                animationBox.max() = max(animationBox.max(), p);
            }
            animationFile<< endl;
        }
    }
    Info << endl;

    // Expand the bound box
    {
        const vector c = animationBox.midpoint(), d = animationBox.span();
        if (d.x() < d.y()/2)
        {
            animationBox.min().x() = c.x() - d.y()/4;
            animationBox.max().x() = c.x() + d.y()/4;
        }
        else if (d.y() < d.x()/2)
        {
            animationBox.min().y() = c.y() - d.x()/4;
            animationBox.max().y() = c.y() + d.x()/4;
        }
    }

    // Write the animation commands
    animationFile
        << "end" << endl << endl
        << "set size ratio -1" << endl
        << "do for [i=1:" << label(endTime/deltaT - 0.5) << "] {" << endl
        << "    set title sprintf('\%g s', i*" << deltaT << ')' << endl
        << "    plot ["
        << animationBox.min().x() << ':' << animationBox.max().x() << "]["
        << animationBox.min().y() << ':' << animationBox.max().y() << ']';
    forAll(outlines, outlineI)
    {
        const label bodyID = outlineBodyIDs[outlineI];
        animationFile<< (outlineI ? "," : "") << " \\" << endl
            << "        '$data' us 2:3 every :::"
            << outlines.size() << "*i+" << outlineI << "::"
            << outlines.size() << "*i+" << outlineI << " w l lw 2 "
            << "t '" << motion.bodies()[bodyID].name() << '\'';
    }
    animationFile<< endl << "    pause " << deltaT << endl << "}" << endl;

    return 0;
}


// ************************************************************************* //
