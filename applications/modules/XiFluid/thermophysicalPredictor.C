/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

    const volScalarField fuDot
    (
        bSource/(b + 0.001)
      - 2*Db*c
       *mag(fvc::grad(ft))/max(ft, scalar(1e-6))
       *mag(fvc::grad(fuFres))/max(fuFres, scalar(1e-6))
    );

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

    // Burnt gas density
    const volScalarField rhob("rhob", thermo.rhob());

    // Calculate flame normal etc.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volVectorField n("n", fvc::grad(b));

    volScalarField mgb("mgb", mag(n));

    const dimensionedScalar dMgb
    (
        "dMgb",
        mgbCoeff_*
        (b*c*mgb)().weightedAverage(mesh.V())
       /((b*c)().weightedAverage(mesh.V()) + small)
      + dimensionedScalar(mgb.dimensions(), small)
    );

    // Stabilise mgb for division here and sub-models
    mgb = max(mgb, dMgb);

    const surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());
    surfaceVectorField nfVec(fvc::interpolate(n));
    nfVec += SfHat*(fvc::snGrad(b) - (SfHat & nfVec));
    nfVec /= max(mag(nfVec), dMgb);
    const surfaceScalarField nf("nf", mesh.Sf() & nfVec);
    n /= mgb;

    // Ignition corrected Xi*Su
    const volScalarField SuXiCorr(Su*this->XiCorr(Xi, nf, dMgb));

    // Turbulent flame speed volumetric flux
    const surfaceScalarField phivSt("phivSt", fvc::interpolate(SuXiCorr)*nf);

    // Turbulent flame speed mass flux
    const surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*SuXiCorr)*nf);

    const dimensionedScalar mgbStab(2*dMgb/bMin_/mgbCoeff_);

    fvScalarMatrix bSourceEqn
    (
        fvModels().source(rho, b)
      - fvm::div(phiSt, b)
      + fvm::Sp(fvc::div(phiSt), b)
      - fvm::Sp(rhou*Su*Xi*mgbStab*max(bMin_ - b, scalar(0)), b)
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

    // Unburnt-mean and burnt-mean gas velocity flux differences
    // caused by the heat-release dilatation.
    //
    // The bounded non-conservative part of the unburnt-mean gas velocity flux
    // difference is already included in the b-Xi burning rate but it can also
    // be included in the Hau equation in bounded non-conserative form for
    // consistency
    //
    // This term cannot be used in partially-premixed combustion due to
    // the resultant inconsistency between ft and Hau transport.
    // A possible solution would be to solve for fuu as well as ft.
    tmp<surfaceScalarField> phiDbu;
    tmp<surfaceScalarField> phib;
    if (!thermo_.containsSpecie("ft"))
    {
        const surfaceScalarField bf(fvc::interpolate(b));

        const surfaceScalarField phiDb
        (
            fvc::interpolate(Db)*fvc::snGrad(b)*mesh.magSf()
        );

        phiDbu = new surfaceScalarField("phiDbu", - phiDb/max(bf, 0.0001));

        // Not currently required
        // phiDbb = new surfaceScalarField("phiDbb", phiDb/max(1 - bf, 0.0001));
    }

    // Correct the flame wrinkling
    XiModel_->correct();

    // Correct the laminar flame speed
    SuModel_->correct();

    if (thermo_.containsSpecie("fu"))
    {
        fuSolve(mvConvection, Db, (bSourceEqn & b));
    }

    HauSolve(mvConvection, Db, phiDbu);
}


void Foam::solvers::XiFluid::HauSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db,
    const tmp<surfaceScalarField>& phiDbu
)
{
    volScalarField& hau = thermo_.heu();

    const volScalarField::Internal rhoByRhou(rho()/thermo.rhou()());

    fvScalarMatrix HauEqn
    (
        fvm::ddt(rho, hau) + mvConvection.fvmDiv(phi, hau)
      + rhoByRhou
       *(
            (fvc::ddt(rho, K) + fvc::div(phi, K))()
          + pressureWork(-dpdt)
        )
      - fvm::laplacian(Db, hau)
     ==
        fvModels().source(rho, hau)
    );

    if (phiDbu.valid())
    {
        HauEqn +=
            fvm::div(phiDbu(), hau, "div(phiSt,b)")
          - fvm::Sp(fvc::div(phiDbu()), hau);
    }

    HauEqn.relax();
    fvConstraints().constrain(HauEqn);
    HauEqn.solve();
    fvConstraints().constrain(hau);
}


void Foam::solvers::XiFluid::HaSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& ha = thermo_.he();

    fvScalarMatrix HaEqn
    (
        fvm::ddt(rho, ha) + mvConvection.fvmDiv(phi, ha)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + pressureWork(-dpdt)
      - fvm::laplacian(Db, ha)
     ==
        (
            buoyancy.valid()
          ? fvModels().source(rho, ha) + rho*(U & buoyancy->g)
          : fvModels().source(rho, ha)
        )
    );

    HaEqn.relax();
    fvConstraints().constrain(HaEqn);
    HaEqn.solve();
    fvConstraints().constrain(ha);
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

    HaSolve(mvConvection(), Db);

    if (!ignited)
    {
        thermo_.heu() == thermo.he();
    }

    thermo_.correct();
}


// ************************************************************************* //
