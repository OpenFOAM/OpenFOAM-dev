/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "argList.H"

#include "vector.H"
#include "IFstream.H"

#include "BSpline.H"
#include "CatmullRomSpline.H"

using namespace Foam;

inline Ostream& printPoint(Ostream& os, const point& p)
{
    os  << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("file .. fileN");
    argList::addBoolOption("B", "B-Spline implementation");
    argList::addBoolOption("CMR", "catmull-rom spline (default)");
    argList::addOption
    (
        "n",
        "INT",
        "number of segments for evaluation - default 20"
    );

    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        args.printUsage();
    }

    bool useBSpline = args.optionFound("B");
    bool useCatmullRom = args.optionFound("CMR");
    label nSeg = args.optionLookupOrDefault<label>("n", 20);

    if (!useCatmullRom && !useBSpline)
    {
        Info<<"defaulting to Catmull-Rom spline" << endl;
        useCatmullRom = true;
    }

    for (label argI=1; argI < args.size(); ++argI)
    {
        const string& srcFile = args[argI];
        Info<< nl << "reading " << srcFile << nl;
        IFstream ifs(srcFile);

        List<pointField> pointFields(ifs);


        forAll(pointFields, splineI)
        {
            Info<<"\n# points:" << endl;
            forAll(pointFields[splineI], ptI)
            {
                printPoint(Info, pointFields[splineI][ptI]);
            }

            if (useBSpline)
            {
                BSpline spl(pointFields[splineI]);

                Info<< nl << "# B-Spline" << endl;

                for (label segI = 0; segI <= nSeg; ++segI)
                {
                    scalar lambda = scalar(segI)/scalar(nSeg);
                    printPoint(Info, spl.position(lambda));
                }
            }

            if (useCatmullRom)
            {
                CatmullRomSpline spl(pointFields[splineI]);

                Info<< nl <<"# Catmull-Rom" << endl;

                for (label segI = 0; segI <= nSeg; ++segI)
                {
                    scalar lambda = scalar(segI)/scalar(nSeg);
                    printPoint(Info, spl.position(lambda));
                }
            }

            {
                polyLine pl(pointFields[splineI]);

                Info<< nl <<"# polyList" << endl;

                for (label segI = 0; segI <= nSeg; ++segI)
                {
                    scalar lambda = scalar(segI)/scalar(nSeg);
                    printPoint(Info, pl.position(lambda));
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
