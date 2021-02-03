/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Test-Function1

Description
    Tests Function1

\*---------------------------------------------------------------------------*/

#include "Function1.H"
#include "IFstream.H"
#include "OFstream.H"
#include "ListOps.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("dictionary");
    argList args(argc, argv);
    const word dictName(args[1]);
    Info<< "Reading " << dictName << nl << endl;
    const dictionary dict = IFstream(dictName)();

    Info<< "Constructing the function\n" << endl;
    const autoPtr<Function1<scalar>> functionPtr =
        Function1<scalar>::New("function", dict);
    const Function1<scalar>& function = functionPtr();

    const scalar x0 = dict.lookup<scalar>("x0");
    const scalar x1 = dict.lookup<scalar>("x1");
    const label nX = dict.lookup<label>("nX");
    const scalar dx = (x1 - x0)/(nX - 1);
    const scalarField xs
    (
        x0 + dx*List<scalar>(identity(nX))
    );

    Info<< "Calculating values\n" << endl;
    const scalarField ys(function.value(xs));
    const scalarField integralYs(function.integral(scalarField(nX, x0), xs));
    scalarField trapezoidIntegralYs(nX, 0);
    for (label i = 1; i < nX; ++ i)
    {
        trapezoidIntegralYs[i] =
            trapezoidIntegralYs[i-1] + (ys[i] + ys[i-1])/2*dx;
    }

    OFstream dataFile(dictName + ".dat");
    Info<< "Writing to " << dataFile.name() << "\n" << endl;
    dataFile<< "# x y integralY trapezoidIntegralY" << endl;
    forAll(xs, i)
    {
        dataFile
            << xs[i] << ' '
            << ys[i] << ' '
            << integralYs[i] << ' '
            << trapezoidIntegralYs[i] << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
