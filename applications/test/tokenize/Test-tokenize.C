/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Test the tokenizing of various things
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "cpuTime.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("string .. stringN");
    argList::addOption("file", "name");
    argList::addOption("repeat", "count");

    argList args(argc, argv, false, true);

    const label repeat = args.optionLookupOrDefault<label>("repeat", 1);

    cpuTime timer;
    for (label count = 0; count < repeat; ++count)
    {
        for (label argI=1; argI < args.size(); ++argI)
        {
            const string& rawArg = args[argI];
            if (count == 0)
            {
                Info<< "input string: " << rawArg << nl;
            }

            IStringStream is(rawArg);

            while (is.good())
            {
                token tok(is);
                // char ch;
                // is.get(ch);
                // is.putback(ch);
                int lookahead = is.peek();

                if (count == 0)
                {
                    Info<< "token: " << tok.info();
                    Info<< "  lookahead: '" << char(lookahead) << "'" << endl;
                }
            }

            if (count == 0)
            {
                Info<< nl;
                IOobject::writeDivider(Info);
            }
        }
    }

    Info<< "tokenized args " << repeat << " times in "
        << timer.cpuTimeIncrement() << " s\n\n";

    if (args.optionFound("file"))
    {
        for (label count = 0; count < repeat; ++count)
        {
            IFstream is(args["file"]);

            if (count == 0)
            {
                Info<< "tokenizing file: " << args["file"] << nl;
            }

            while (is.good())
            {
                token tok(is);
                if (count == 0)
                {
                    Info<< "token: " << tok.info() << endl;
                }
            }

            if (count == 0)
            {
                Info<< nl;
                IOobject::writeDivider(Info);
            }
        }

        Info<< "tokenized file " << repeat << " times in "
            << timer.cpuTimeIncrement() << " s\n\n";
    }

    return 0;
}

// ************************************************************************* //
