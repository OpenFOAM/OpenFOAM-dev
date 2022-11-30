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
    applyBoundaryLayer

Description
    Apply a simplified boundary-layer model to the velocity and
    turbulence fields based on the 1/7th power-law.

    The uniform boundary-layer thickness is either provided via the -ybl option
    or calculated as the average of the distance to the wall scaled with
    the thickness coefficient supplied via the option -Cbl.  If both options
    are provided -ybl is used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "wallDist.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Turbulence constants - file-scope
static const scalar Cmu(0.09);
static const scalar kappa(0.41);


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "apply a simplified boundary-layer model to the velocity and\n"
        "turbulence fields based on the 1/7th power-law."
    );

    argList::addOption
    (
        "ybl",
        "scalar",
        "specify the boundary-layer thickness"
    );
    argList::addOption
    (
        "Cbl",
        "scalar",
        "boundary-layer thickness as Cbl * mean distance to wall"
    );
    argList::addBoolOption
    (
        "writenut",
        "write nut field"
    );

    #include "setRootCase.H"

    if (!args.optionFound("ybl") && !args.optionFound("Cbl"))
    {
        FatalErrorInFunction
            << "Neither option 'ybl' or 'Cbl' have been provided to calculate "
            << "the boundary-layer thickness.\n"
            << "Please choose either 'ybl' OR 'Cbl'."
            << exit(FatalError);
    }
    else if (args.optionFound("ybl") && args.optionFound("Cbl"))
    {
        FatalErrorInFunction
            << "Both 'ybl' and 'Cbl' have been provided to calculate "
            << "the boundary-layer thickness.\n"
            << "Please choose either 'ybl' OR 'Cbl'."
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createMeshNoChangers.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Modify velocity by applying a 1/7th power law boundary-layer
    // u/U0 = (y/ybl)^(1/7)
    // assumes U0 is the same as the current cell velocity

    Info<< "Setting boundary layer velocity" << nl << endl;
    scalar yblv = ybl.value();
    forAll(U, celli)
    {
        if (y[celli] <= yblv)
        {
            mask[celli] = 1;
            U[celli] *= ::pow(y[celli]/yblv, (1.0/7.0));
        }
    }
    mask.correctBoundaryConditions();

    Info<< "Writing U\n" << endl;
    U.write();

    // Update/re-write phi
    #include "createPhi.H"
    phi.write();

    autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

    autoPtr<incompressible::momentumTransportModel> turbulence
    (
        incompressible::momentumTransportModel::New(U, phi, viscosity)
    );

    if (isA<incompressible::RASModel>(turbulence()))
    {
        // Calculate nut
        turbulence->validate();
        tmp<volScalarField> tnut = turbulence->nut();
        volScalarField& nut = const_cast<volScalarField&>(tnut());
        volScalarField S(mag(dev(symm(fvc::grad(U)))));
        nut = (1 - mask)*nut + mask*sqr(kappa*min(y, ybl))*::sqrt(2)*S;

        // Do not correct BC - wall functions will 'undo' manipulation above
        // by using nut from turbulence model

        if (args.optionFound("writenut"))
        {
            Info<< "Writing nut" << endl;
            nut.write();
        }


        //--- Read and modify turbulence fields

        // Turbulence k
        tmp<volScalarField> tk = turbulence->k();
        volScalarField& k = const_cast<volScalarField&>(tk());
        scalar ck0 = pow025(Cmu)*kappa;
        k = (1 - mask)*k + mask*sqr(nut/(ck0*min(y, ybl)));
        k.correctBoundaryConditions();

        Info<< "Writing k\n" << endl;
        k.write();


        // Turbulence epsilon
        tmp<volScalarField> tepsilon = turbulence->epsilon();
        volScalarField& epsilon = const_cast<volScalarField&>(tepsilon());
        scalar ce0 = ::pow(Cmu, 0.75)/kappa;
        epsilon = (1 - mask)*epsilon + mask*ce0*k*sqrt(k)/min(y, ybl);

        // Do not correct BC - G set by the wall-functions is not available
        // epsilon.correctBoundaryConditions();

        Info<< "Writing epsilon\n" << endl;
        epsilon.write();

        // Turbulence omega
        typeIOobject<volScalarField> omegaHeader
        (
            "omega",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (omegaHeader.headerOk())
        {
            volScalarField omega(omegaHeader, mesh);

            const incompressible::RASModel& rasModel =
                refCast<const incompressible::RASModel>(turbulence());

            omega = (1 - mask)*omega + mask*ce0*sqrt(k)/(Cmu*min(y, ybl));
            bound(omega, rasModel.omegaMin());

            // Do not correct BC - G set by the wall-functions is not available
            // omega.correctBoundaryConditions();

            Info<< "Writing omega\n" << endl;
            omega.write();
        }

        // Turbulence nuTilda
        typeIOobject<volScalarField> nuTildaHeader
        (
            "nuTilda",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (nuTildaHeader.headerOk())
        {
            volScalarField nuTilda(nuTildaHeader, mesh);
            nuTilda = nut;

            // Do not correct BC - G set by the wall-functions is not available
            // nuTilda.correctBoundaryConditions();

            Info<< "Writing nuTilda\n" << endl;
            nuTilda.write();
        }
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
