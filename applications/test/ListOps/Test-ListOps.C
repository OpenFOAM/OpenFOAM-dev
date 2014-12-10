/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Test-ListOps

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "SubList.H"
#include "ListOps.H"
#include "face.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Test Rotations:" << nl << endl;

    List<label> forwardRotate(identity(5));
    face testFace(identity(4));

    for (label i = 0; i < 8; ++i)
    {
        Info<< "Rotate forward by " << i << " : "
            << rotateList(forwardRotate, i) << endl;
    }

    for (label i = 0; i < 8; ++i)
    {
        Info<< "Rotate backward by " << i << " : "
            << rotateList(forwardRotate, -i) << endl;
    }

    Info<< nl << "Face                 : " << testFace << endl;
    Info<< "Rotate by 2          : " << rotateList(testFace, 2) << endl;
    inplaceRotateList<List, label>(testFace, -6);
    Info<< "Rotate inplace by -6 : " << testFace << nl << endl;

    Info<< "Test inplace rotate      : " << forwardRotate << endl;
    inplaceRotateList(forwardRotate, 2);
    Info<< "Rotate to the right by 2 : " << forwardRotate << endl;
    inplaceRotateList(forwardRotate, -2);
    Info<< "Rotate to the left by 2  : " << forwardRotate << endl;

    List<label> subRotate(identity(10));
    SubList<label> subL(subRotate, 5, 3);

    Info<< "Test inplace rotate on sublist : " << subRotate << endl;
    inplaceRotateList(subL, 3);
    Info<< "Rotate to the right by 3       : " << subRotate << endl;
    inplaceRotateList(subL, -8);
    Info<< "Rotate to the left by 3        : " << subRotate << endl;

    Info<< nl << nl << "Test Reversing:" << nl << endl;

    Info<< "List    : " << identity(5) << endl;
    Info<< "Reverse : " << reverseList(identity(5)) << endl;
    Info<< "List    : " << identity(6) << endl;
    Info<< "Reverse : " << reverseList(identity(6)) << nl << endl;

    List<label> test1(identity(5));
    Info<< "List            : " << test1 << endl;
    inplaceReverseList(test1);
    Info<< "Inplace Reverse : " << test1 << nl << endl;

    List<label> test2(identity(6));
    Info<< "List            : " << test2 << endl;
    inplaceReverseList(test2);
    Info<< "Inplace Reverse : " << test2 << nl << endl;

    face test3(identity(6));
    Info<< "Face            : " << test3 << endl;
    inplaceReverseList(test3);
    Info<< "Inplace Reverse : " << test3 << nl << endl;

    FixedList<label, 6> test4(identity(6));
    Info<< "FixedList       : " << test4 << endl;
    inplaceReverseList(test4);
    Info<< "Inplace Reverse : " << test4 << nl << endl;

    List<label> test5(identity(9));
    SubList<label> test5SubList(test5, 4, 3);
    Info<< "List                            : " << test5 << endl;
    inplaceReverseList(test5SubList);
    Info<< "Reverse Sublist between 3 and 6 : " << test5 << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
