/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Test-quaternion

Description
    Test application for quaternions.

\*---------------------------------------------------------------------------*/

#include "quaternion.H"
#include "septernion.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    quaternion q(vector(1, 0, 0), 0.7853981);
    Info<< "q " << q << endl;

    vector v(0, 1, 0);
    Info<< "v " << v << endl;

    Info<< "inv(q)*q " << inv(q)*q << endl;

    Info<< "q*quaternion(0, v)*conjugate(q) "
        << q*quaternion(0, v)*conjugate(q) << endl;

    Info<< "q.transform(v) " << q.transform(v) << endl;
    Info<< "q.R() & v " << (q.R() & v) << endl;

    Info<< "q.invTransform(v) " << q.invTransform(v) << endl;

    septernion tr(vector(0, 0.1, 0), q);
    Info<< "tr " << tr << endl;

    Info<< "inv(tr)*tr " << inv(tr)*tr << endl;

    Info<< "tr.transform(v) " << tr.transform(v) << endl;

    Info<< "(septernion(vector(0, -1, 0))*q*septernion(vector(0, 1, 0)))"
        << ".transform(v) "
        <<  (septernion(vector(0, -1, 0))
           *q
           *septernion(vector(0, 1, 0))).transform(v)
        << endl;

    Info<< "Test conversion from and to Euler-angles" << endl;
    vector angles(0.1, 0.2, 0.3);
    for (int rs=quaternion::ZYX; rs<quaternion::XZX; ++rs)
    {
        if
        (
            mag
            (
                angles
              - quaternion(quaternion::rotationSequence(rs), angles)
               .eulerAngles(quaternion::rotationSequence(rs))
            )
          > SMALL
        )
        {
            FatalErrorInFunction
                << "Inconsistent conversion for rotation sequence "
                << rs << exit(FatalError)
                << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
