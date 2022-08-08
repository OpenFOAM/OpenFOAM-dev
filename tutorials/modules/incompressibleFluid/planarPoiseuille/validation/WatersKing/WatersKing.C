/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    WatersKing

Description
    Analytical solution for the start-up planar Poiseuille flow of an
    Oldroyd-B fluid.

    References:
    \verbatim
        Waters, N. D., & King, M. J. (1970).
        Unsteady flow of an elasto-viscous liquid.
        Rheologica Acta, 9, 345-355.

        Amoreira, L. J., & Oliveira, P. J. (2010).
        Comparison of different formulations for the numerical
        calculation of unsteady incompressible viscoelastic fluid
        flow. Adv. Appl. Math. Mech, 4, 483-502.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "createMesh.H"
    #include "createFields.H"

    const scalar h = mesh.bounds().span().y();
    Info<< "Height from centreline to wall = " << h << endl;

    label centrelineID = mesh.boundary().findPatchID("centreline");
    const vector patchToCell =
        mesh.boundary()[centrelineID].Cf()[0]
      - mesh.C()[mesh.findNearestCell(location)];

    const scalar y = patchToCell.y()/h;
    Info<< "Normalised distance from centreline = " << y << nl << endl;

    const scalar nu0 = nu1 + nu2;
    const scalar E = lambda*nu0/(rho*sqr(h));
    const scalar beta = nu2/nu0;
    const scalar UInf = K*sqr(h)/3.0/nu0;

    Info<< "Waters and King parameters:" << nl
        << "E =    " << E << nl
        << "beta = " << beta << nl
        << "K =    " << K << nl
        << "UInf = " << UInf << nl << endl;

    label order = 8;

    scalarField ak(order, 0);
    scalarField bk(order, 0);
    scalarField ck(order, 0);
    scalarField B(order, 0);

    forAll(ak, i)
    {
        scalar k = i + 1;
        ak[i] = (2.0*k - 1)/2.0*constant::mathematical::pi*::sqrt(E);
        bk[i] = (1.0 + beta*sqr(ak[i]))/2.0;
        ck[i] = ::sqrt(mag(sqr(bk[i]) - sqr(ak[i])));
        B[i]  = 48*::pow(-1, k)/::pow((2*k - 1)*constant::mathematical::pi, 3)*
                ::cos((2*k - 1)*constant::mathematical::pi*y/2);
    }

    scalarField A(order, 0);
    OFstream file(runTime.path()/"WatersKing.dat");
    const scalar LOGvGreat = ::log(vGreat);
    while (!runTime.end())
    {
        scalar t = runTime.userTimeValue()/lambda;
        forAll(A, i)
        {
            if (bk[i]*t < LOGvGreat)
            {
                if (bk[i] >= ak[i])
                {
                    A[i] = (bk[i] - sqr(ak[i]))/ck[i]*::sinh(ck[i]*t)
                    + ::cosh(ck[i]*t);
                }
                else
                {
                    A[i] = (bk[i] - sqr(ak[i]))/ck[i]*::sin(ck[i]*t)
                         + ::cos(ck[i]*t);
                }
                A[i] *= ::exp(-bk[i]*t);
            }
            else
            {
                Info<< "Coefficient A[" << order << "] = 0" << endl;
                order = i;
                Info<< "Resizing A and B to " << order << endl;
                A.resize(order);
                B.resize(order);
            }
        }
        scalar U = UInf*(1.5*(1 - sqr(y)) + sum(A*B));
        file<< runTime.timeName() << token::TAB << U << endl;
        runTime++;
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
