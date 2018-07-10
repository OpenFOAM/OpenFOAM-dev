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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "IStringStream.H"
#include "labelList.H"
#include "DynamicList.H"
#include "face.H"

using namespace Foam;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    List<label> lstA(10);
    List<label> lstC(IStringStream("(1 2 3 4)")());

    forAll(lstA, i)
    {
        lstA[i] = i;
    }

    Info<< "lstA: " << lstA << endl;
    Info<< "lstC: " << lstC << endl;

    Xfer<List<label>> xA = xferMove(lstA);
    Xfer<List<label>> xB;

    List<label> lstB( xA );

    Info<< "xA: " << xA() << endl;
    Info<< "xB: " << xB() << endl;
    Info<< "lstA: " << lstA << endl;
    Info<< "lstB: " << lstB << endl;
    Info<< "lstC: " << lstC << endl;

    xA = lstB;

    Info<< "xA: " << xA() << endl;
    Info<< "xB: " << xB() << endl;
    Info<< "lstA: " << lstA << endl;
    Info<< "lstB: " << lstB << endl;
    Info<< "lstC: " << lstC << endl;

    xB = xA;

    List<label> lstD(xferCopy(lstC));
    List<label> lstE(xferMove(lstC));

    // this must be empty
    List<label> lstF = xferCopy(lstC);

    Info<< "xA: " << xA() << endl;
    Info<< "xB: " << xB() << endl;
    Info<< "lstA: " << lstA << endl;
    Info<< "lstB: " << lstB << endl;
    Info<< "lstC: " << lstC << endl;
    Info<< "lstD: " << lstD << endl;
    Info<< "lstE: " << lstE << endl;
    Info<< "lstF: " << lstF << endl;

    Info<< "xB[" << xB->size() << "]\n";

    // clear the underlying List
    xB->clear();

    Info<< "xB[" << xB->size() << "]\n";

    DynamicList<label> dl(10);
    for (label i = 0; i < 5; ++i)
    {
        dl.append(i);
    }

    face f1(dl);
    face f2(xferCopy<labelList>(dl));

    Info<< "dl[" << dl.size() << "/" << dl.capacity() << "] " << dl << endl;
    Info<< "f1: " << f1 << endl;
    Info<< "f2: " << f2 << endl;

    // add some more labels
    for (label i = 5; i < 8; ++i)
    {
        dl.append(i);
    }

    // note: xfer() method returns a plain labelList
    face f3(dl.xfer());
    Info<< "dl[" << dl.size() << "/" << dl.capacity() << "] " << dl << endl;
    Info<< "f3: " << f3 << endl;

    Info<<"\nflip faces:" << endl;
    f1.flip();
    f3.flip();
    Info<< "f1: " << f1 << endl;
    Info<< "f3: " << f3 << endl;

    return 0;
}


// ************************************************************************* //
