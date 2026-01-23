/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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
#include "fvcFlux.H"

#include "EulerDdtScheme.H"
#include "gaussConvectionScheme.H"
#include "upwind.H"

#include "CMULES.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::XiFluid::burn()
{
    volScalarField& b(thermo_.b());

    // Progress variable
    volScalarField& c(thermo_.c());

    // Unburnt gas density
    const volScalarField& rhou(thermo_.uThermo().rho());

    const volScalarField Db(XiModel_->Db());

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

    // Turbulent flame speed volumetric flux
    const surfaceScalarField phivSt("phivSt", fvc::interpolate(Su*Xi)*nf);

    // Turbulent flame speed mass flux
    surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*Su*Xi)*nf);

    const dimensionedScalar mgbStab(2*dMgb/bMin_/mgbCoeff_);

    surfaceScalarField phib("phib", phi + phiSt);

    tmp<surfaceScalarField> tbPhiUD;
    tmp<surfaceScalarField> tbPhiStUD;
    tmp<surfaceScalarField> tbLaplacianPhi;
    tmp<surfaceScalarField> tbLaplacianPhiCorr;
    tmp<volScalarField::Internal> tSu;
    tmp<volScalarField::Internal> tSp;

    // Bounded implicit b predictor
    {
        // Construct the bounded upwind interpolator for b
        const upwind<scalar> bUD(mesh, phib);

        // Construct the corrected Laplacian
        fvScalarMatrix bLaplacian(fvm::laplacian(Db, b));

        // Set the Laplacian correction source to 0
        bLaplacian.source() = Zero;

        // Construct the b convection matrix
        fvScalarMatrix bPhi
        (
            fv::gaussConvectionScheme<scalar>(mesh, phib, bUD).fvmDiv(phi, b)
        );

        // Construct the b flame propagation matrix
        fvScalarMatrix bPhiSt
        (
            fv::gaussConvectionScheme<scalar>(mesh, phib, bUD).fvmDiv(phiSt, b)
        );

        const volScalarField::Internal divPhiSt(fvc::div(phiSt));

        //- Construct the b source matrix
        fvScalarMatrix bSource
        (
            fvModels().source(rho, b)
          - fvm::Sp(rhou*Su*Xi*mgbStab*max(bMin_ - b, scalar(0)), b)
        );

        // Assemble the bounded b matrix
        fvScalarMatrix bEqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(rho, b)
          + bPhi + bPhiSt
          - bLaplacian
         ==
            bSource
          + fvm::Sp(divPhiSt, b)
        );

        // Solve for b and constrain
        bEqn.relax();
        fvConstraints().constrain(bEqn);
        bEqn.solve();
        fvConstraints().constrain(b);

        // Set the fluxes for the MULES corrector

        tbPhiUD = bPhi.flux();
        tbPhiStUD = bPhiSt.flux();

        tbLaplacianPhi = -bLaplacian.flux();
        if (bLaplacian.faceFluxCorrectionPtr())
        {
            tbLaplacianPhiCorr = -*bLaplacian.faceFluxCorrectionPtr();
        }

        tSu = bSource.Su() + divPhiSt*b;
        tSp = bSource.Sp();
    }

    const volScalarField::Internal& Sp = tSp();

    const word divbName("div(phi,b)");

    // Higher-order face interpolate of b
    const surfaceScalarField bf
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            phib,
            mesh.schemes().div(divbName)
        )().interpolate(phib, b)
    );

    // Higher-order convection and flame-propagation face-fluxes
    const surfaceScalarField bPhi(phi*bf);
    const surfaceScalarField bPhiSt(phiSt*bf);

    // Higher-order convection and flame-propagation face-flux corrections
    surfaceScalarField bPhiCorr(bPhi - tbPhiUD());
    surfaceScalarField bPhiStCorr(bPhiSt - tbPhiStUD());

    if (tbLaplacianPhiCorr.valid())
    {
        bPhiCorr += tbLaplacianPhiCorr();
    }

    const MULES::control MULEScontrols(mesh.solution().solverDict(b.name()));

    // Cache a list of the flux corrections to be limited
    UPtrList<surfaceScalarField> bPhiCorrs{&bPhiCorr, &bPhiStCorr};

    // MULES limited bounded explicit b corrector
    MULES::correct
    (
        MULEScontrols,
        rho,
        b,
        tbPhiUD() + tbPhiStUD(),
        bPhiCorrs,
        Sp,
        oneField(),
        zeroField()
    );

    // Recalculate c from b
    c = scalar(1) - b;

    // Correct the flame wrinkling
    XiModel_->correct();

    // Correct the laminar flame speed
    SuModel_->correct();

    // Set the b-equation convection+diffusion flux
    // for the solution of the unburnt gas energy and species
    phib = tbPhiUD() + tbLaplacianPhi() + bPhiCorr;

    // Set the c-equation convection+diffusion flux
    // for the solution of the burnt gas energy and species
    const surfaceScalarField phic("phic", phi - phib);

    // Set the b-equation source term
    // for the solution of the unburnt and burnt gas energy and species
    const volScalarField::Internal bSource
    (
        tSu() + tSp()*b()
      - fvc::div(tbPhiStUD() + bPhiStCorr)()
    );

    // Set the unburnt and burnt gas stabilisation coefficients
    const volScalarField::Internal bStab(max(bMin_ - b, scalar(0)));
    const volScalarField::Internal cStab(max(bMin_ - c, scalar(0)));

    // Solve for the unburnt gas energy and species
    uSolve(bStab, phib, bSource);

    // Solve for the burnt gas energy and species
    bSolve(cStab, phic, bSource);
}


void Foam::solvers::XiFluid::uSolve
(
    const volScalarField::Internal& bStab,
    const surfaceScalarField& phib,
    const volScalarField::Internal& bSource
)
{
    PtrList<volScalarField>& Yu = thermo_.uThermo().Y();

    uReaction_->correct();

    forAll(Yu, i)
    {
        volScalarField& Yui = Yu[i];

        if (uThermo.solveSpecie(i))
        {
            uSolve
            (
                Yui,
                Yu.size() > 2 ? "Yi" : Yui.name(),
                bStab,
                phib,
                fvm::Sp(bSource, Yui)
            );
        }
        else
        {
            Yui.correctBoundaryConditions();
        }
    }

    thermo_.uThermo().normaliseY();

    HuSolve(bStab, phib, bSource);
}


void Foam::solvers::XiFluid::bSolve
(
    const volScalarField::Internal& cStab,
    const surfaceScalarField& phic,
    const volScalarField::Internal& bSource
)
{
    PtrList<volScalarField>& Yb = thermo_.bThermo().Y();

    if (Yb.size())
    {
        bReaction_->correct();

        const PtrList<volScalarField::Internal> Yp(uThermo.prompt());

        forAll(Yb, i)
        {
            volScalarField& Ybi = Yb[i];

            if (bThermo.solveSpecie(i))
            {
                bSolve
                (
                    Ybi,
                    Yb.size() > 2 ? "Yi" : Ybi.name(),
                    cStab,
                    phic,
                    fvm::Su(-bSource*Yp[i], Ybi)
                );
            }
            else
            {
                Ybi.correctBoundaryConditions();
            }
        }

        thermo_.bThermo().normaliseY();
    }

    HbSolve(cStab, phic, bSource);
}


Foam::tmp<Foam::fvScalarMatrix> Foam::solvers::XiFluid::fvmStab
(
    const volScalarField& bc,
    const volScalarField::Internal& bcStab,
    const volScalarField& D,
    volScalarField& f
)
{
    // Advective-diffusive stabilisation for b,c -> 0
    return bcStab*
    (
        fv::EulerDdtScheme<scalar>(mesh).fvmDdt(rho, f)
      + fv::gaussConvectionScheme<scalar>
        (
            mesh,
            phi,
            upwind<scalar>(mesh, phi)
        ).fvmDiv(phi, f)
      - fv::gaussLaplacianScheme<scalar, scalar>::fvmLaplacianUncorrected
        (
            fvc::interpolate(D)*mesh.magSf(),
            mesh.nonOrthDeltaCoeffs(),
            f
        )
    );
}


void Foam::solvers::XiFluid::ubSolve
(
    volScalarField& f,
    const word& fName,
    const volScalarField& alpha,
    const volScalarField& bc,
    const volScalarField::Internal& bcStab,
    const surfaceScalarField& phibc,
    const volScalarField& D,
    const thermophysicalTransportModel& thermophysicalTransport,
    const fvScalarMatrix& combustionRate,
    const reactionModel& reaction
)
{
    fvScalarMatrix fEqn
    (
        fvm::ddt(bc, rho, f)
      + fvm::div(phibc, f, "div(" + phibc.name() + ',' + fName + ')')

        // Advective-diffusive stabilisation for bc -> 0
      + fvmStab(bc, bcStab, D, f)

        // Diffusive transport within the unburnt/burnt gas
      + thermophysicalTransport.divj(f)
     ==
        // Combustion source
        combustionRate

        // Reaction rate within the unburnt/burnt gas
      + alpha*reaction.R(f)

        // Other sources
      + fvModels().source(bc, rho, f)
    );

    fEqn.relax();
    fvConstraints().constrain(fEqn);
    fEqn.solve(fName);
    fvConstraints().constrain(f);
}


void Foam::solvers::XiFluid::uSolve
(
    volScalarField& fu,
    const word& fuName,
    const volScalarField::Internal& bStab,
    const surfaceScalarField& phib,
    const fvScalarMatrix& source
)
{
    const volScalarField Du("Du", rho*(momentumTransport.nut() + uThermo.nu()));
    ubSolve
    (
        fu,
        fuName,
        thermo_.alphau(),
        b,
        bStab,
        phib,
        Du,
        uThermophysicalTransport_(),
        source,
        uReaction_()
    );
}


void Foam::solvers::XiFluid::bSolve
(
    volScalarField& fb,
    const word& fbName,
    const volScalarField::Internal& cStab,
    const surfaceScalarField& phic,
    const fvScalarMatrix& source
)
{
    const volScalarField Db("Db", rho*(momentumTransport.nut() + bThermo.nu()));
    ubSolve
    (
        fb,
        fbName,
        thermo_.alphab(),
        c,
        cStab,
        phic,
        Db,
        bThermophysicalTransport_(),
        source,
        bReaction_()
    );
}


void Foam::solvers::XiFluid::HuSolve
(
    const volScalarField::Internal& bStab,
    const surfaceScalarField& phib,
    const volScalarField::Internal& bSource
)
{
    volScalarField& hu = thermo_.uThermo().he();

    const volScalarField::Internal rhoByRhou(rho()/uThermo.rho()());

    const volScalarField Du("Du", rho*(momentumTransport.nut() + uThermo.nu()));

    fvScalarMatrix HuEqn
    (
        fvm::ddt(b, rho, hu) + fvm::div(phib, hu)

        // Advective-diffusive stabilisation for b -> 0
      + fvmStab(b, bStab, Du, hu)

        // Pressure-work
      + fvc::ddt(b, rho, K) + fvc::div(phib, K)
      + (b + bStab)*rhoByRhou*pressureWork(-dpdt)

        // Diffusive transport within the unburnt gas
      + uThermophysicalTransport_->divq(hu)
     ==
        // Combustion source
        fvm::Sp(bSource, hu)

        // Other sources
      + fvModels().source(b, rho, hu)
    );

    HuEqn.relax();
    fvConstraints().constrain(HuEqn);
    HuEqn.solve();
    fvConstraints().constrain(hu);
}


void Foam::solvers::XiFluid::HbSolve
(
    const volScalarField::Internal& cStab,
    const surfaceScalarField& phic,
    const volScalarField::Internal& bSource
)
{
    volScalarField& hb = thermo_.bThermo().he();

    const volScalarField::Internal rhoByRhob(rho()/bThermo.rho()());

    const volScalarField Db("Db", rho*(momentumTransport.nut() + bThermo.nu()));

    fvScalarMatrix HbEqn
    (
        fvm::ddt(c, rho, hb) + fvm::div(phic, hb)

        // Advective-diffusive stabilisation for b -> 0
      + fvmStab(c, cStab, Db, hb)

        // Pressure-work
      + fvc::ddt(c, rho, K) + fvc::div(phic, K)
      + (c + cStab)*rhoByRhob*pressureWork(-dpdt)

        // Diffusive transport within the unburnt gas
      + uThermophysicalTransport_->divq(hb)
     ==
        // Combustion source
      - bSource*(uThermo.ha()() - bThermo.hf()())()

        // Other sources
      + fvModels().source(c, rho, hb)
    );

    HbEqn.relax();
    fvConstraints().constrain(HbEqn);
    HbEqn.solve();
    fvConstraints().constrain(hb);
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

    // At the point of ignition initialise the burnt gas
    // thermophysical properties
    if (ignited && !ignited_)
    {
        ignited_ = ignited;

        if (thermo_.bThermo().Y().size())
        {
            const PtrList<volScalarField::Internal> Yp(uThermo.prompt());

            // Approximate phic for Ybi boundary condition correction
            const surfaceScalarField phic("phic", phi);

            forAll(Yp, i)
            {
                thermo_.bThermo().Y(i).internalFieldRef() = Yp[i];
                thermo_.bThermo().Y(i).correctBoundaryConditions();
            }
        }

        thermo_.bThermo().he() = uThermo.ha() - bThermo.hf();
        thermo_.bThermo().correct();
    }

    if (ignited_)
    {
        burn();
    }
    else
    {
        const volScalarField::Internal bStab
        (
            IOobject
            (
                "bStab",
                runTime.time().name(),
                mesh
            ),
            mesh,
            scalar(0)
        );

        const surfaceScalarField phib("phib", phi);

        const volScalarField::Internal bSource
        (
            IOobject
            (
                "bSource",
                runTime.time().name(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimDensity/dimTime, 0)
        );

        uSolve(bStab, phib, bSource);
    }

    thermo_.correct();
}


// ************************************************************************* //
