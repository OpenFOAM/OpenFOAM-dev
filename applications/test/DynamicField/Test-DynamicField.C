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

\*---------------------------------------------------------------------------*/

#include "point.H"
#include "DynamicField.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    DynamicField<point, 0, 10, 11> testField;
    DynamicField<point, 0, 10, 11> testField2;

    testField.setSize(5);
    testField2.setSize(5);

    testField[0] = testField2[0] = vector(1.0, 4.5, 6.3);
    testField[1] = testField2[1] = vector(5.2, 2.3, 3.5);
    testField[2] = testField2[2] = vector(7.5, 4.7, 7.7);
    testField[3] = testField2[3] = vector(2.8, 8.2, 2.3);
    testField[4] = testField2[4] = vector(6.1, 1.7, 8.8);

    Info << "testField:" << testField << endl;

    testField.append(vector(0.5, 4.8, 6.2));

    Info << "testField after appending:" << testField << endl;

    testField.append(vector(2.7, 2.3, 6.1));

    Info << "testField after appending:" << testField << endl;

    vector elem = testField.remove();

    Info << "removed element:" << elem << endl;
    Info << "testField:" << testField << endl;

    testField.append(vector(3.0, 1.3, 9.2));

    Info << "testField:" << testField << endl;

    testField.setSize(10, vector(1.5, 0.6, -1.0));

    Info << "testField after setSize:" << testField << endl;

    testField.append(testField2);

    Info << "testField after appending testField2:" << testField << endl;

    testField = testField2;

    Info << "testField after assignment:" << testField << endl;

    testField += testField2;

    Info << "testField after field algebra:" << testField << endl;

    testField.clear();

    testField.append(vector(3.0, 1.3, 9.2));

    Info << "testField after clear and append:" << testField << endl;

    testField.clearStorage();

    Info << "testField after clearStorage:" << testField << endl;

    return 0;
}


// ************************************************************************* //
