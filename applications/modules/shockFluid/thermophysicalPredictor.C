/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "shockFluid.H"
#include "fvmDdt.H"
#include "fvcDiv.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockFluid::thermophysicalPredictor()
{
    volScalarField& e = thermo_.he();

    const surfaceScalarField e_pos(interpolate(e, pos, thermo.T().name()));
    const surfaceScalarField e_neg(interpolate(e, neg, thermo.T().name()));

    surfaceScalarField phiEp
    (
        "phiEp",
        aphiv_pos()*(rho_pos()*(e_pos + 0.5*magSqr(U_pos())) + p_pos())
      + aphiv_neg()*(rho_neg()*(e_neg + 0.5*magSqr(U_neg())) + p_neg())
      + aSf()*(p_pos() - p_neg())
    );

    // Make flux for pressure-work absolute
    if (mesh.moving())
    {
        phiEp += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
    }

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvc::div(phiEp)
      + fvc::ddt(rho, K)
     ==
        fvModels().source(rho, e)
    );

    if (!inviscid)
    {
        const surfaceScalarField devTauDotU
        (
            "devTauDotU",
            devTau() & (a_pos()*U_pos() + a_neg()*U_neg())
        );

        EEqn += thermophysicalTransport->divq(e) + fvc::div(devTauDotU);
    }

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(e);

    thermo_.correct();
}


// ************************************************************************* //
