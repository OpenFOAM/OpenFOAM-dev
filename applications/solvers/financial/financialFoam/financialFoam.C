/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    financialFoam

Description
    Solves the Black-Scholes equation to price commodities.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include "setWriter.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating value(price of commodities)" << endl;

    surfaceScalarField phi("phi", (sigmaSqr - r)*(Pf & mesh.Sf()));

    volScalarField DV("DV", 0.5*sigmaSqr*sqr(P.component(Foam::vector::X)));

    Info<< "Starting time loop\n" << endl;

    while (runTime.loop())
    {
        delta == fvc::grad(V)().component(Foam::vector::X);

        solve
        (
            fvm::ddt(V)
          + fvm::div(phi, V)
          - fvm::Sp(fvc::div(phi), V)
          - fvm::laplacian(DV, V)
         ==
          - fvm::Sp(r, V)
        );

        runTime.write();

        if (runTime.writeTime())
        {
            setWriter::New(runTime.graphFormat())->write
            (
                runTime.globalPath()
               /functionObjects::writeFile::outputPrefix
               /args.executable()
               /runTime.name(),

                args.executable(),

                coordSet(true, word::null, mesh.C().primitiveField(), "x"),

                "V", V.primitiveField(),
                "delta", delta.primitiveField()
            );
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
