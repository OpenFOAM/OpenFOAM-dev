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
    Test some string functionality

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "stringOps.H"
#include "dictionary.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    string test
    (
        "  $HOME kjhkjhkjh \" \\$HOME/tyetyery $; ${FOAM_RUN} \n $; hkjh;"
        " $(DONOTSUBST) some other <${USER}> with '${__UNKNOWN:-some default}'"
        " value "
        " or with '${HOME:+Home was set}' via :+ alternative"
        " or with '${__UNKNOWN:+unknown}' empty"
    );

    dictionary dict;
    dict.add("HOME", "myHome");

    dictionary subDict;
    subDict.add("value1", "test1");
    subDict.add("value2", "test2");
    dict.add("FOAM_RUN", subDict);


    Info<< "string:" << test << nl << "hash:"
        << unsigned(string::hash()(test)) << endl;

    Info<<"trimLeft: " << stringOps::trimLeft(test) << endl;
    Info<<"trimRight: " << stringOps::trimRight(test) << endl;
    Info<<"trim: " << stringOps::trim(test) << endl;


    // test sub-strings via iterators
    string::const_iterator iter  = test.end();
    string::const_iterator iter2 = test.end();
    string::size_type fnd = test.find('\\');

    if (fnd != string::npos)
    {
        iter  = test.begin() + fnd;
        iter2 = iter + 6;
    }

    Info<< "sub-string via iterators : >";
    while (iter != iter2)
    {
        Info<< *iter;
        ++iter;
    }
    Info<< "<\n";

    Info<< string(test).replace("kj", "zzz") << endl;
    Info<< string(test).replace("kj", "") << endl;
    Info<< string(test).replaceAll("kj", "zzz") << endl;
    Info<< string(test).replaceAll("kj", "z") << endl;

    Info<< "expanded: " << string(test).expand() << endl;

    Info<<"dictionary-based substitution: " << dict << endl;
    Info<< "expand dict: " << stringOps::expand(test, dict) << endl;

    string test2("~OpenFOAM/controlDict");
    Info<< test2 << " => " << test2.expand() << endl;

    // replace controlDict with "newName"
    {
        string::size_type i = test2.rfind('/');

        if (i == string::npos)
        {
            test2 = "newName";
        }
        else
        {
            // this uses the std::string::replace
            test2.replace(i+1, string::npos, "newName");
        }
        Info<< "after replace: " << test2 << endl;

        // do another replace
        // this uses the Foam::string::replace
        test2.replace("OpenFOAM", "openfoam");

        Info<< "after replace: " << test2 << endl;
    }

    string s;
    Sin.getLine(s);

    string s2(s.expand());

    cout<< "output string with " << s2.length() << " characters\n";
    cout<< "ostream<<  >" << s2 << "<\n";
    Info<< "Ostream<<  >" << s2 << "<\n";
    Info<< "hash:" << hex << string::hash()(s2) << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
