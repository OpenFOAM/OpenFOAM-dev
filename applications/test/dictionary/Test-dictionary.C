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
    dictionaryTest

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "stringOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("dict .. dictN");
    argList args(argc, argv, false, true);

    {
        dictionary dict;
        dict.add("aaOPENMPIcc", 1);

        string s("DDD${aa${WM_MPLIB}cc}EEE");
        stringOps::inplaceExpand(s, dict, true, false);
        Info<< "variable expansion:" << s << endl;
    }


    Info<< nl
        << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
        << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
        << endl;

    if (args.size() <= 1)
    {
        {
            dictionary dict1(IFstream("testDict")());
            Info<< "dict1: " << dict1 << nl
                << "toc: " << dict1.toc() << nl
                << "keys: " << dict1.keys() << nl
                << "patterns: " << dict1.keys(true) << endl;

            dictionary dict2(dict1.xfer());

            Info<< "dict1.toc(): " << dict1.name() << " " << dict1.toc() << nl
                << "dict2.toc(): " << dict2.name() << " " << dict2.toc()
                << endl;

            // copy back
            dict1 = dict2;
            Info<< "dict1.toc(): " << dict1.name() << " " << dict1.toc()
                << endl;

            dictionary dict3(dict2.subDictPtr("boundaryField"));
            dictionary dict4(dict2.subDictPtr("NONEXISTENT"));

            Info<< "dictionary construct from pointer" << nl
                << "ok = " << dict3.name() << " " << dict3.toc() << nl
                << "no = " << dict4.name() << " " << dict4.toc() << endl;
        }


        IOobject::writeDivider(Info);

        {
            dictionary dict(IFstream("testDictRegex")());
            dict.add(keyType("fooba[rz]", true), "anything");

            Info<< "dict:" << dict << nl
                << "toc: " << dict.toc() << nl
                << "keys: " << dict.keys() << nl
                << "patterns: " << dict.keys(true) << endl;

            Info<< "Pattern find \"abc\" in top directory : "
                << dict.lookup("abc") << endl;
            Info<< "Pattern find \"abc\" in sub directory : "
                << dict.subDict("someDict").lookup("abc")
                    << endl;
            Info<< "Recursive pattern find \"def\" in sub directory : "
                << dict.subDict("someDict").lookup("def", true)
                    << endl;
            Info<< "Recursive pattern find \"foo\" in sub directory : "
                << dict.subDict("someDict").lookup("foo", true)
                    << endl;
            Info<< "Recursive pattern find \"fooz\" in sub directory : "
                << dict.subDict("someDict").lookup("fooz", true)
                    << endl;
            Info<< "Recursive pattern find \"bar\" in sub directory : "
                << dict.subDict("someDict").lookup("bar", true)
                    << endl;
            Info<< "Recursive pattern find \"xxx\" in sub directory : "
                << dict.subDict("someDict").lookup("xxx", true)
                    << endl;
        }
    }
    else
    {
        IOobject::writeDivider(Info);
        for (label argI=1; argI < args.size(); ++argI)
        {
            const string& dictFile = args[argI];
            IFstream is(dictFile);

            dictionary dict(is);

            Info<< dict << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
