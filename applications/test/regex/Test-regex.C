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

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "regExp.H"
#include "List.H"
#include "Tuple2.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    List<Tuple2<string, string> > rawList(IFstream("testRegexps")());
    Info<< "input list:" << rawList << endl;
    IOobject::writeDivider(Info) << endl;

    List<string> groups;

    // report matches:
    forAll(rawList, elemI)
    {
        const string& pat = rawList[elemI].first();
        const string& str = rawList[elemI].second();
        regExp re(pat);

        Info<< str << " =~ m/" << pat.c_str() << "/ == ";

        if (re.match(str, groups))
        {
            Info<< "true";
            if (re.ngroups())
            {
                Info<< " groups:" << groups;
            }
        }
        else
        {
            Info<< "false";
            if (re.search(str))
            {
                Info<< " partial match";
            }
        }
        Info<< endl;
    }

    Info<<"test regExp(const char*) ..." << endl;
    string me("Mark");

    if (regExp("[Mm]ar[ck]").match(me))
    {
        Info<< "matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    if (regExp("").match(me))
    {
        Info<< "matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    if (regExp(NULL).match(me))
    {
        Info<< "matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    Info<< endl;

    return 0;
}


// ************************************************************************* //
