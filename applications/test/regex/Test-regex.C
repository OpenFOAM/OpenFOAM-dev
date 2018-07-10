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

Description
    Tests for regular expressions

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
    List<Tuple2<string, string>> rawList(IFstream("testRegexps")());
    Info<< "Test expressions:" << rawList << endl;
    IOobject::writeDivider(Info) << endl;

    List<string> groups;

    // Report matches:
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
                Info<< nl << "groups: " << groups;
            }
        }
        else
        {
            if (re.search(str))
            {
                Info<< " partial match";
            }
            else
            {
                Info<< "false";
            }
        }
        Info<< endl;
    }

    Info<< nl << "test regExp(const char*) ..." << endl;
    string me("Mark");

    // Handling of null strings
    if (regExp(nullptr).match(me))
    {
        Info<< "fail - matched: " << me << endl;
    }
    else
    {
        Info<< "pass - null pointer is no expression" << endl;
    }

    // Normal match
    if (regExp("[Mm]ar[ck]").match(me))
    {
        Info<< "pass - matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    // Match ignore case
    if (regExp("mar[ck]", true).match(me))
    {
        Info<< "pass - matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    // Embedded prefix for match ignore case
    if (regExp("(?i)mar[ck]").match(me))
    {
        Info<< "pass - matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    // Handling of empty expression
    if (regExp("").match(me))
    {
        Info<< "fail - matched: " << me << endl;
    }
    else
    {
        Info<< "pass - no match on empty expression" << endl;
    }

    // Embedded prefix - but expression is empty
    if (regExp("(?i)").match(me))
    {
        Info<< "fail - matched: " << me << endl;
    }
    else
    {
        Info<< "pass - no match on empty expression" << endl;
    }


    Info<< endl;

    return 0;
}


// ************************************************************************* //
