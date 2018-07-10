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

#include "dimensionedTensor.H"
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    dimensionedTensor dt("dt", dimLength, tensor(0, 1, 2, 3, 4, 5, 6, 7, 8));

    Info<< "dt.component(tensor::XX): " << dt.component(tensor::XX) << endl;

    dimensionedScalar ds("ds", dimTime, 1.0);

    Info<< "ds*dt dt*ds: " << ds*dt << " " << dt*ds << endl;

    dimensionedTensor dt2("dt2", dimLength, tensor(1, 1, 2, 3, 4, 5, 6, 7, 8));

    Info<< "cmptMultiply(dt, dt2): " << cmptMultiply(dt, dt2) << endl;
    Info<< "cmptDivide(dt, dt2): " << cmptDivide(dt, dt2) << endl;

    {
        Pout<< "dimensionSet construct from is:"
            << dimensionSet(IStringStream("[Pa m^2 s^-2]")())
            << endl;

        IStringStream is("[Pa m^2 s^-2]");
        dimensionSet dset(dimless);
        is >> dset;
        Pout<< "dimensionSet read:" << dset << endl;
    }

    {
        Pout<< "construct from is:"
            << dimensionedScalar(IStringStream("bla [Pa mm^2 s^-2] 3.0")())
            << endl;
        Pout<< "construct from name,is:"
            <<  dimensionedScalar
                (
                    "ABC",
                    IStringStream("[Pa mm^2 s^-2] 3.0")()
                ) << endl;
        Pout<< "construct from name,dimensionSet,is:"
            <<  dimensionedScalar
                (
                    "ABC",
                    dimLength,
                    IStringStream("bla [mm] 3.0")()
                ) << endl;
        {
            IStringStream is("bla [mm] 3.0");
            dimensionedScalar ds;
            is >> ds;
            Pout<< "read:" << ds << endl;
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
