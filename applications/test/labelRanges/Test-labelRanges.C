/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Description
    Test label ranges
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "labelRanges.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("start size .. startN sizeN");
    argList::addOption("verbose");
    argList::addNote
    (
        "The default is to add ranges, use 'add' and 'del' to toggle\n\n"
        "Eg, 0 10 30 10 del 20 15"
    );

    argList args(argc, argv, false, true);

    if (args.optionFound("verbose"))
    {
        labelRange::debug = 1;
    }


    labelRanges ranges;

    bool removeMode = false;
    for (label argI=1; argI < args.size()-1; ++argI)
    {
        if (args[argI] == "add")
        {
            removeMode = false;
            continue;
        }
        else if (args[argI] == "del")
        {
            removeMode = true;
            continue;
        }

        label start = 0;
        label size  = 0;

        IStringStream(args[argI])() >> start;
        ++argI;
        IStringStream(args[argI])() >> size;

        labelRange range(start, size);

        Info<< "---------------" << nl;
        if (removeMode)
        {
            Info<< "del " << range << " :";
            forAllConstIter(labelRange, range, iter)
            {
                Info<< " " << iter();
            }
            Info<< nl;

            ranges.remove(range);
        }
        else
        {
            Info<< "add " << range  << " :";
            forAllConstIter(labelRange, range, iter)
            {
                Info<< " " << iter();
            }
            Info<< nl;

            ranges.add(range);
        }

        Info<< "<list>" << ranges << "</list>" << nl;
        forAllConstIter(labelRanges, ranges, iter)
        {
            Info<< " " << iter();
        }
        Info<< nl;
    }

    return 0;
}

// ************************************************************************* //
