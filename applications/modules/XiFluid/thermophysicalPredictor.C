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
#include "fvcDdt.H"
#include "fvmDiv.H"

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

    fvConstraints().constrain(ftEqn);

    ftEqn.solve();

    fvConstraints().constrain(ft);
}


Foam::tmp<Foam::volScalarField> Foam::solvers::XiFluid::XiCorr
(
    const volScalarField& c,
    const surfaceScalarField& nf,
    const dimensionedScalar& dMgb,
    const volScalarField& Xi
) const
{
    // Calculate volume of ignition kernel
    const dimensionedScalar Vk
    (
        "Vk",
        dimVolume,
        gSum(c*mesh.V().primitiveField())
    );

    // Radius of the ignition kernel
    dimensionedScalar rk("rk", dimLength, 0.0);

    // Area of the ignition kernel
    dimensionedScalar Ak("Ak", dimArea, 0.0);

    if (Vk.value() > small)
    {
        // Calculate kernel area from its volume
        // and the dimensionality of the case

        switch(mesh.nGeometricD())
        {
            case 3:
            {
                // Assume it is part-spherical
                const dimensionedScalar sphereFraction
                (
                    "ignitionSphereFraction",
                    dimless,
                    combustionProperties
                );

                rk = pow
                    (
                        (3.0/4.0)*Vk
                       /(sphereFraction*constant::mathematical::pi),
                        1.0/3.0
                    );

                Ak = sphereFraction*4*constant::mathematical::pi*sqr(rk);
            }
            break;

            case 2:
            {
                // Assume it is part-cylindrical
                const dimensionedScalar thickness
                (
                    "ignitionThickness",
                    dimLength,
                    combustionProperties
                );

                const dimensionedScalar circleFraction
                (
                    "ignitionCircleFraction",
                    dimless,
                    combustionProperties
                );

                rk = sqrt
                (
                    Vk
                   /(
                       circleFraction
                      *constant::mathematical::pi
                      *thickness
                    )
                );

                Ak = circleFraction*2*constant::mathematical::pi*rk
                    *thickness;
            }
            break;

            case 1:
                // Assume it is planar with given area
                Ak = dimensionedScalar
                (
                    "ignitionKernelArea",
                    dimArea,
                    combustionProperties
                );

                rk = Vk/Ak;
            break;
        }

        const dimensionedScalar maxXiCorrRadius
        (
            "maxKernelCorrRadius",
            dimLength,
            combustionProperties
        );

        if (rk.value() < maxXiCorrRadius.value())
        {
            // Calculate kernel area from b field consistent with the
            // discretisation of the b equation.
            const volScalarField mgb
            (
                fvc::div(nf, b, "div(phiSt,b)") - b*fvc::div(nf) + dMgb
            );
            const dimensionedScalar AkEst
            (
                gSum(mgb*mesh.V().primitiveField())
            );

            const scalar XiCorr = max(min((Ak/AkEst).value(), 10.0), 1.0);
            Info<< "XiCorr = " << XiCorr << endl;

            return XiCorr*Xi;
        }
    }

    return Xi;
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
        fvc::interpolate(rhou*Su*XiCorr(c, nf, dMgb, Xi))*nf
    );

    // Create b equation
    fvScalarMatrix bEqn
    (
        fvm::ddt(rho, b)
      + mvConvection.fvmDiv(phi, b)
      + fvm::div(phiSt, b)
      - fvm::Sp(fvc::div(phiSt), b)
      - fvm::laplacian(Db, b)
     ==
        fvModels().source(rho, b)
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
}


void Foam::solvers::XiFluid::EauSolve
(
    const fv::convectionScheme<scalar>& mvConvection,
    const volScalarField& Db
)
{
    volScalarField& heau = thermo_.heu();

    const volScalarField::Internal rhoByRhou(rho()/thermo.rhou()());

    fvScalarMatrix heauEqn
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

    fvConstraints().constrain(heauEqn);

    heauEqn.solve();

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
    ignited_ =
        mesh().time().value() - mesh().time().deltaTValue() >= ignitionStart_;

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
        ignited_ ? XiModel_->Db() : thermophysicalTransport.DEff(b)
    );

    if (thermo_.containsSpecie("ft"))
    {
        ftSolve(mvConvection(), Db);
    }

    if (ignited_)
    {
        bSolve(mvConvection(), Db);
        EauSolve(mvConvection(), Db);
    }

    EaSolve(mvConvection(), Db);

    if (!ignited_)
    {
        thermo_.heu() == thermo.he();
    }

    thermo_.correct();
}


// ************************************************************************* //
