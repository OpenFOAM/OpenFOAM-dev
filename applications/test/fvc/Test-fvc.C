/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    test

Description
    Finite volume method test code.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    volScalarField fx(pow(mesh.C().component(vector::X), 1));
    fx.write();
    volScalarField gradx4(fvc::grad(fx)().component(vector::X));
    gradx4.write();

    volVectorField curlC(fvc::curl(1.0*mesh.C()));
    curlC.write();

    surfaceScalarField xf(mesh.Cf().component(vector::X));
    surfaceScalarField xf4(pow(xf, 4));

    for (int i=1; i<xf4.size()-1; i++)
    {
        scalar gradx4a = (xf4[i] - xf4[i-1])/(xf[i] - xf[i-1]);
        Info<< (gradx4a - gradx4[i])/gradx4a << endl;
    }

    Info<< "end" << endl;
}


// ************************************************************************* //
