/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "XiFluid.H"
#include "bXiIgnition.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvcSup.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::XiFluid::ftSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& ft = thermo_.Y("ft");

    fvScalarMatrix ftEqn
    (
        fvm::ddt(rho, ft)
      + mvConvection.fvmDiv(phi, ft)
      - fvm::laplacian(Db, ft)
     ==
        fvModels().source(rho, ft)
    );

    ftEqn.relax();
    fvConstraints().constrain(ftEqn);
    ftEqn.solve();
    fvConstraints().constrain(ft);
}


void Foam::solvers::XiFluid::fuSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db,
    const volScalarField& bSource
)
{
    volScalarField& fu = thermo_.Y("fu");
    const volScalarField& ft = thermo_.Y("ft");
    const volScalarField& b(b_);

    // Progress variable
    const volScalarField c("c", scalar(1) - b);

    // Unburnt gas density
    const volScalarField rhou("rhou", thermo.rhou());

    const volScalarField fres(thermo_.fres());
    const volScalarField fuFres(max(fu - fres, scalar(0)));

    const volScalarField fuDot =
        bSource/(b + 0.001)
      - 2*Db*c
       *mag(fvc::grad(ft))/max(ft, 1e-6)
       *mag(fvc::grad(fuFres))/max(fuFres, 1e-6);

    fvScalarMatrix fuEqn
    (
        fvm::ddt(rho, fu)
      + mvConvection.fvmDiv(phi, fu)
      - fvm::laplacian(Db, fu)
     ==
        fvm::Sp(fuDot, fu)
      - fuDot*fres
      + fvModels().source(rho, fu)
    );

    fuEqn.relax();
    fvConstraints().constrain(fuEqn);
    fuEqn.solve();
    fvConstraints().constrain(fu);
}


void Foam::solvers::XiFluid::egrSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& egr = thermo_.Y("egr");

    fvScalarMatrix egrEqn
    (
        fvm::ddt(rho, egr)
      + mvConvection.fvmDiv(phi, egr)
      - fvm::laplacian(Db, egr)
     ==
        fvModels().source(rho, egr)
    );

    egrEqn.relax();
    fvConstraints().constrain(egrEqn);
    egrEqn.solve();
    fvConstraints().constrain(egr);
}



Foam::tmp<Foam::volScalarField> Foam::solvers::XiFluid::XiCorr
(
    const volScalarField& Xi,
    const surfaceScalarField& nf,
    const dimensionedScalar& dMgb
) const
{
    const UPtrListDictionary<fv::bXiIgnition> ignitionModels
    (
        fvModels().lookupType<fv::bXiIgnition>()
    );

    bool igniting = false;

    forAll(ignitionModels, i)
    {
        if (ignitionModels[i].igniting())
        {
            igniting = true;
            break;
        }
    }

    if (igniting)
    {
        tmp<volScalarField> tXi(volScalarField::New("XiCorrected", Xi));

        // Calculate kernel area from b field consistent with the
        // discretisation of the b equation.
        const volScalarField mgb
        (
            fvc::div(nf, b, "div(phiSt,b)") - b*fvc::div(nf) + dMgb
        );

        forAll(ignitionModels, i)
        {
            ignitionModels[i].XiCorr(tXi.ref(), b, mgb);
        }

        return tXi;
    }
    else
    {
        return Xi;
    }
}


void Foam::solvers::XiFluid::bSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& b(b_);

    // Progress variable
    const volScalarField c("c", scalar(1) - b);

    // Unburnt gas density
    const volScalarField rhou("rhou", thermo.rhou());

    // Calculate flame normal etc.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volVectorField n("n", fvc::grad(b));

    volScalarField mgb("mgb", mag(n));

    const dimensionedScalar dMgb
    (
        "dMgb",
        1.0e-3*
        (b*c*mgb)().weightedAverage(mesh.V())
       /((b*c)().weightedAverage(mesh.V()) + small)
      + dimensionedScalar(mgb.dimensions(), small)
    );

    mgb += dMgb;

    const surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());
    surfaceVectorField nfVec(fvc::interpolate(n));
    nfVec += SfHat*(fvc::snGrad(b) - (SfHat & nfVec));
    nfVec /= (mag(nfVec) + dMgb);
    surfaceScalarField nf("nf", mesh.Sf() & nfVec);
    n /= mgb;

    // Turbulent flame speed flux
    const surfaceScalarField phiSt
    (
        "phiSt",
        fvc::interpolate(rhou*Su*XiCorr(Xi, nf, dMgb))*nf
    );

    fvScalarMatrix bSourceEqn
    (
        fvModels().source(rho, b)
      - fvm::div(phiSt, b)
      + fvm::Sp(fvc::div(phiSt), b)
    );

    // Create b equation
    fvScalarMatrix bEqn
    (
        fvm::ddt(rho, b)
      + mvConvection.fvmDiv(phi, b)
      - fvm::laplacian(Db, b)
     ==
        bSourceEqn
    );

    // Solve for b and constrain
    bEqn.relax();
    fvConstraints().constrain(bEqn);
    bEqn.solve();
    fvConstraints().constrain(b);

    // Correct the flame wrinkling
    XiModel_->correct();

    // Correct the laminar flame speed
    SuModel_->correct();

    if (thermo_.containsSpecie("fu"))
    {
        fuSolve(mvConvection, Db, (bSourceEqn & b));
    }
}


void Foam::solvers::XiFluid::EauSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& heau = thermo_.heu();

    const volScalarField::Internal rhoByRhou(rho()/thermo.rhou()());

    fvScalarMatrix EauEqn
    (
        fvm::ddt(rho, heau) + mvConvection.fvmDiv(phi, heau)
      + rhoByRhou
       *(
            (fvc::ddt(rho, K) + fvc::div(phi, K))()
          + pressureWork
            (
                heau.name() == "eau"
              ? mvConvection.fvcDiv(phi, p/rho)()
              : -dpdt
            )
        )
      - fvm::laplacian(Db, heau)

        // These terms cannot be used in partially-premixed combustion due to
        // the resultant inconsistency between ft and heau transport.
        // A possible solution would be to solve for ftu as well as ft.
        //- fvm::div(muEff*fvc::grad(b)/(b + 0.001), heau)
        //+ fvm::Sp(fvc::div(muEff*fvc::grad(b)/(b + 0.001)), heau)

     ==
        fvModels().source(rho, heau)
    );

    EauEqn.relax();
    fvConstraints().constrain(EauEqn);
    EauEqn.solve();
    fvConstraints().constrain(heau);
}


void Foam::solvers::XiFluid::EaSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& hea = thermo_.he();

    fvScalarMatrix EaEqn
    (
        fvm::ddt(rho, hea) + mvConvection.fvmDiv(phi, hea)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + pressureWork
        (
            hea.name() == "ea"
          ? mvConvection.fvcDiv(phi, p/rho)()
          : -dpdt
        )
      - fvm::laplacian(Db, hea)
     ==
        (
            buoyancy.valid()
          ? fvModels().source(rho, hea) + rho*(U & buoyancy->g)
          : fvModels().source(rho, hea)
        )
    );

    EaEqn.relax();
    fvConstraints().constrain(EaEqn);
    EaEqn.solve();
    fvConstraints().constrain(hea);
}


void Foam::solvers::XiFluid::thermophysicalPredictor()
{
    const UPtrListDictionary<fv::bXiIgnition> ignitionModels
    (
        fvModels().lookupType<fv::bXiIgnition>()
    );

    bool ignited = false;

    forAll(ignitionModels, i)
    {
        if (ignitionModels[i].ignited())
        {
            ignited = true;
            break;
        }
    }

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.schemes().div("div(phi,ft_b_ha_hau)")
        )
    );

    const volScalarField Db
    (
        "Db",
        ignited ? XiModel_->Db() : thermophysicalTransport.DEff(b)
    );

    if (thermo_.containsSpecie("ft"))
    {
        ftSolve(mvConvection(), Db);
    }

    if (ignited)
    {
        bSolve(mvConvection(), Db);
        EauSolve(mvConvection(), Db);
    }
    else
    {
        if (thermo_.containsSpecie("fu"))
        {
            volScalarField& fu = thermo_.Y("fu");
            const volScalarField& ft = thermo_.Y("ft");
            fu = ft;
        }
    }

    if (thermo_.containsSpecie("egr"))
    {
        egrSolve(mvConvection(), Db);
    }

    EaSolve(mvConvection(), Db);

    if (!ignited)
    {
        thermo_.heu() == thermo.he();
    }

    thermo_.correct();
}


// ************************************************************************* //
