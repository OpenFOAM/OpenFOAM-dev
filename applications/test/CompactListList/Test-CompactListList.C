/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    CompactListListTest

Description
    Simple demonstration and test application for the CompactListList class.

\*---------------------------------------------------------------------------*/

#include "CompactListList.H"
#include "IOstreams.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "faceList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    {
        // null construct
        CompactListList<label> cll1;
        Info<< "cll1:" << cll1 << endl;

        // Resize and assign row by row
        labelList row0(2, 0);
        labelList row1(3, 1);

        labelList rowSizes(2);
        rowSizes[0] = row0.size();
        rowSizes[1] = row1.size();
        cll1.resize(rowSizes);

        cll1[0].deepCopy(row0);
        cll1[1].deepCopy(row1);
        Info<< "cll1:" << cll1 << endl;

        forAll(cll1.m(), i)
        {
            Info<< "i:" << i << " whichRow:" << cll1.whichRow(i) << endl;
        }
    }

    List<List<label>> lll(5);
    lll[0].setSize(3, 0);
    lll[1].setSize(2, 1);
    lll[2].setSize(6, 2);
    lll[3].setSize(0, 3);
    lll[4].setSize(1, 4);

    CompactListList<label> cll2(lll);

    Info<< "cll2  = " << cll2 << endl;

    forAll(cll2, i)
    {
        Info<< cll2[i] << endl;
    }

    Info<< endl;

    Info<< "cll2(2, 3) = " << cll2(2, 3) << nl << endl;
    cll2(2, 3) = 999;
    Info<< "cll2(2, 3) = " << cll2(2, 3) << nl << endl;

    Info<< "cll2 as List<List<label >> " << cll2()
        << endl;

    cll2.setSize(3);

    Info<< "cll2  = " << cll2 << endl;

    cll2.setSize(0);

    Info<< "cll2  = " << cll2 << endl;


    List<label> rowSizes(5);
    rowSizes[0] = 2;
    rowSizes[1] = 0;
    rowSizes[2] = 1;
    rowSizes[3] = 3;
    rowSizes[4] = 2;

    CompactListList<label> cll3(rowSizes, 1);

    Info<< "cll3 = " << cll3 << endl;

    CompactListList<label> cll4;

    cll4.transfer(cll3);

    Info<< "cll3 = " << cll3 << endl;
    Info<< "cll4 = " << cll4 << endl;


    {
        // IO
        OStringStream ostr;
        ostr << cll4;

        IStringStream istr(ostr.str());
        CompactListList<label> cll5(istr);
        Info<< "cll5 = " << cll5 << endl;
    }
    {
        // IO
        cll4.clear();
        OStringStream ostr;
        ostr << cll4;

        IStringStream istr(ostr.str());
        CompactListList<label> cll5(istr);
        Info<< "cll5 = " << cll5 << endl;
    }

    {
        faceList fcs(2);
        fcs[0] = face(labelList(1, 111));
        fcs[1] = face(labelList(2, 222));

        CompactListList<label, face> compactFcs(fcs);
        Info<< "compactFcs:" << compactFcs << endl;

        faceList fcs2 = compactFcs();
        Info<< "fcs2:" << fcs2 << endl;
    }

    return 0;
}


// ************************************************************************* //
