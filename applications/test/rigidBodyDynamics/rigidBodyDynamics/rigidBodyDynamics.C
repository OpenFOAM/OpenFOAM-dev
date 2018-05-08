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
    const scalar endTime(readScalar(dict.lookup("endTime")));
    const dictionary& bodiesDict = dict.subDict("bodies");
    List<vectorField> bodiesOutline(bodiesDict.size());
    const label i0 = motion.bodies().size() - bodiesOutline.size();
    forAll(bodiesOutline, i)
    {
        const dictionary& bodyDict =
            bodiesDict.subDict(motion.bodies()[i0 + i].name());
        bodiesOutline[i] = vectorField(bodyDict.lookup("outline"));
    }

    // Set up motion fields
    scalarField tau(motion.nDoF(), Zero);
    Field<spatialVector> fx(motion.nBodies(), Zero);
    Info << motion.nDoF() << " degrees of freedom" << endl;

    // Initialise plot bound box and file
    boundBox plotBox;
    OFstream plotFile(dictName + ".gnuplot");
    plotFile<< "$data << end" << endl;

    // Run the RBD simulation
    for (scalar t = 0; t <= endTime + 0.5*deltaT; t += deltaT)
    {
        Info().stdStream() << "\33[2KTime = " << t << '\r' << std::flush;

        motion.newTime();
        motion.solve(t + deltaT, deltaT, tau, fx);

        // Write the bodies' outlines at the current transformation
        forAll(bodiesOutline, i)
        {
            const spatialTransform& invX0 = motion.X0(i0 + i).inv();
            forAll(bodiesOutline[i], j)
            {
                const vector p = invX0.transformPoint(bodiesOutline[i][j]);
                plotFile<< t << ' ' << p.x() << ' ' << p.y() << endl;
                plotBox.min() = min(plotBox.min(), p);
                plotBox.max() = max(plotBox.max(), p);
            }
            plotFile<< endl;
        }
    }
    Info << endl;

    // Write the plot commands
    plotFile
        << "end" << endl << endl
        << "set size ratio -1" << endl
        << "do for [i=1:" << label(endTime/deltaT - 0.5) << "] {" << endl
        << "    set title sprintf('\%g s', i*" << deltaT << ')' << endl
        << "    plot "
        << "[" << plotBox.min().x() << ':' << plotBox.max().x() << ']'
        << "[" << plotBox.min().y() << ':' << plotBox.max().y() << ']';
    forAll(bodiesOutline, i)
    {
        plotFile<< (i ? "," : "") << " \\" << endl
            << "        '$data' us 2:3 every :::"
            << bodiesOutline.size() << "*i+" << i << "::"
            << bodiesOutline.size() << "*i+" << i << " w lp lw 2 pt 7 "
            << "t '" << motion.bodies()[i0 + i].name() << '\'';
    }
    plotFile<< endl << "    pause " << deltaT << endl << "}" << endl;

    return 0;
}


// ************************************************************************* //
