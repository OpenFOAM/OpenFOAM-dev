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
    Test-fileName

Description
    Test some basic fileName functionality

\*---------------------------------------------------------------------------*/

#include "fileName.H"
#include "SubList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "etcFiles.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    wordList wrdList(5);
    wrdList[0] = "hello";
    wrdList[1] = "hello1";
    wrdList[2] = "hello2";
    wrdList[3] = "hello3";
    wrdList[4] = "hello4.hmm";

    fileName pathName(wrdList);

    Info<< "pathName = " << pathName << nl
        << "pathName.name()     = >" << pathName.name() << "<\n"
        << "pathName.path()     = "  << pathName.path() << nl
        << "pathName.ext()      = >" << pathName.ext() << "<\n"
        << "pathName.name(true) = >" << pathName.name(true) << "<\n";

    Info<< "pathName.components() = " << pathName.components() << nl
        << "pathName.component(2) = " << pathName.component(2) << nl
        << endl;

    // try with different combination
    // The final one should emit warnings
    for (label start = 0; start <= wrdList.size(); ++start)
    {
        fileName instance, local;
        word name;

        fileName path(SubList<word>(wrdList, wrdList.size()-start, start));
        fileName path2 = "."/path;

        IOobject::fileNameComponents
        (
            path,
            instance,
            local,
            name
        );

        Info<< "IOobject::fileNameComponents for " << path << nl
            << "  instance = " << instance << nl
            << "  local    = " << local << nl
            << "  name     = " << name << endl;

        IOobject::fileNameComponents
        (
            path2,
            instance,
            local,
            name
        );

        Info<< "IOobject::fileNameComponents for " << path2 << nl
            << "  instance = " << instance << nl
            << "  local    = " << local << nl
            << "  name     = " << name << endl;

    }

    // test findEtcFile
    Info<< "\n\nfindEtcFile tests:" << nl
        << " controlDict => " << findEtcFile("controlDict") << nl
        << " badName => " << findEtcFile("badName") << endl;

    Info<< "This should emit a fatal error:" << endl;
    Info<< " badName(die) => " << findEtcFile("badName", true) << nl
        << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
