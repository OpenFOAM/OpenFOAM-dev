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

\*---------------------------------------------------------------------------*/

#include "stringListOps.H"
#include "IStringStream.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    stringList strLst
    (
        IStringStream
        (
            "("
            "\"hello\""
            "\"heello\""
            "\"heeello\""
            "\"bye\""
            "\"bbye\""
            "\"bbbye\""
            "\"okey\""
            "\"okkey\""
            "\"okkkey\""
            ")"
        )()
    );

    wordReList reLst(IStringStream("( okey \"[hy]e+.*\" )")());

    Info<< "stringList " << strLst << nl;

    labelList matches = findStrings(".*ee.*", strLst);

    Info<< "matches found for regexp .*ee.* :" << nl << matches << nl;
    forAll(matches, i)
    {
        Info<< " -> " << strLst[matches[i]] << nl;
    }
    Info<< endl;

    matches = findStrings(reLst, strLst);

    Info<< "matches found for " << reLst << nl << matches << nl;
    forAll(matches, i)
    {
        Info<< " -> " << strLst[matches[i]] << nl;
    }
    Info<< endl;

    stringList subLst = subsetStrings(".*ee.*", strLst);
    Info<< "subset stringList: " << subLst << nl;

    subLst = subsetStrings(reLst, strLst);
    Info<< "subset stringList: " << subLst << nl;

    inplaceSubsetStrings(reLst, strLst);
    Info<< "subsetted stringList: " << strLst << nl;

    inplaceSubsetStrings(".*l.*", strLst);
    Info<< "subsetted stringList: " << strLst << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
