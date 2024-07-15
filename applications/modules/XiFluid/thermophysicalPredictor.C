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


Foam::dimensionedScalar Foam::solvers::XiFluid::StCorr
(
    const volScalarField& c,
    const surfaceScalarField& nf,
    const dimensionedScalar& dMgb
) const
{
    dimensionedScalar StCorr("StCorr", dimless, 1.0);

    if (ign.igniting())
    {
        // Calculate volume of ignition kernel
        const dimensionedScalar Vk
        (
            "Vk",
            dimVolume,
            gSum(c*mesh.V().primitiveField())
        );
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
                    const scalar sphereFraction
                    (
                        combustionProperties.lookup<scalar>
                        (
                            "ignitionSphereFraction"
                        )
                    );

                    Ak = sphereFraction*4.0*constant::mathematical::pi
                       *pow
                        (
                            3.0*Vk
                           /(sphereFraction*4.0*constant::mathematical::pi),
                            2.0/3.0
                        );
                }
                break;

                case 2:
                {
                    // Assume it is part-circular
                    const dimensionedScalar thickness
                    (
                        combustionProperties.lookup("ignitionThickness")
                    );

                    const scalar circleFraction
                    (
                        combustionProperties.lookup<scalar>
                        (
                            "ignitionCircleFraction"
                        )
                    );

                    Ak = circleFraction*constant::mathematical::pi*thickness
                       *sqrt
                        (
                            4.0*Vk
                           /(
                               circleFraction
                              *thickness
                              *constant::mathematical::pi
                            )
                        );
                }
                break;

                case 1:
                    // Assume it is plane or two planes
                    Ak = dimensionedScalar
                    (
                        combustionProperties.lookup("ignitionKernelArea")
                    );
                break;
            }

            // Calculate kernel area from b field consistent with the
            // discretisation of the b equation.
            const volScalarField mgb
            (
                fvc::div(nf, b, "div(phiSt,b)") - b*fvc::div(nf) + dMgb
            );
            const dimensionedScalar AkEst = gSum(mgb*mesh.V().primitiveField());

            StCorr.value() = max(min((Ak/AkEst).value(), 10.0), 1.0);

            Info<< "StCorr = " << StCorr.value() << endl;
        }
    }

    return StCorr;
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
    const volScalarField rhou(thermo.rhou());

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
        fvc::interpolate(rhou*StCorr(c, nf, dMgb)*Su*Xi)*nf
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


    // Add ignition cell contribution to b-equation
    forAll(ign.sites(), i)
    {
        const ignitionSite& ignSite = ign.sites()[i];

        if (ignSite.igniting())
        {
            forAll(ignSite.cells(), icelli)
            {
                label ignCell = ignSite.cells()[icelli];
                Info<< "Igniting cell " << ignCell;

                Info<< " state :"
                    << ' ' << b[ignCell]
                    << ' ' << Xi[ignCell]
                    << ' ' << Su[ignCell]
                    << ' ' << mgb[ignCell]
                    << endl;

                bEqn.diag()[ignSite.cells()[icelli]] +=
                (
                    ignSite.strength()*ignSite.cellVolumes()[icelli]
                   *rhou[ignSite.cells()[icelli]]/ignSite.duration()
                )/(b[ignSite.cells()[icelli]] + 0.001);
            }
        }
    }

    // Solve for and constrain b
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
        ign.igniting() ? XiModel_->Db() : thermophysicalTransport.DEff(b)
    );

    if (thermo_.containsSpecie("ft"))
    {
        ftSolve(mvConvection(), Db);
    }

    if (ign.ignited())
    {
        bSolve(mvConvection(), Db);
        EauSolve(mvConvection(), Db);
    }

    EaSolve(mvConvection(), Db);

    if (!ign.ignited())
    {
        thermo_.heu() == thermo.he();
    }

    thermo_.correct();
}


// ************************************************************************* //
