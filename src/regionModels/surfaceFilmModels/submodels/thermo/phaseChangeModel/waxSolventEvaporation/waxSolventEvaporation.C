/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "waxSolventEvaporation.H"
#include "thermoSingleLayer.H"
#include "liquidThermo.H"
#include "basicSpecieMixture.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmSup.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waxSolventEvaporation, 0);

addToRunTimeSelectionTable
(
    phaseChangeModel,
    waxSolventEvaporation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar waxSolventEvaporation::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waxSolventEvaporation::waxSolventEvaporation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    speciePhaseChange(typeName, film, dict),
    Wwax_
    (
        IOobject
        (
            IOobject::modelName("Wwax", typeName),
            film.regionMesh().time().constant(),
            film.regionMesh()
        ),
        coeffDict_.lookup<scalar>("Wwax")
    ),
    Wsolvent_
    (
        IOobject
        (
            IOobject::modelName("Wsolvent", typeName),
            film.regionMesh().time().constant(),
            film.regionMesh()
        ),
        coeffDict_.lookup<scalar>("Wsolvent")
    ),
    Ysolvent0_
    (
        IOobject
        (
            IOobject::modelName("Ysolvent0", typeName),
            film.regionMesh().time().constant(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    Ysolvent_
    (
        IOobject
        (
            IOobject::modelName("Ysolvent", typeName),
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    ),
    deltaMin_(coeffDict_.lookup<scalar>("deltaMin")),
    L_(coeffDict_.lookup<scalar>("L")),
    TbFactor_(coeffDict_.lookupOrDefault<scalar>("TbFactor", 1.1)),
    YInfZero_(coeffDict_.lookupOrDefault<Switch>("YInfZero", false)),
    activityCoeff_
    (
        Function1<scalar>::New("activityCoeff", coeffDict_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

waxSolventEvaporation::~waxSolventEvaporation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class YInfType>
void waxSolventEvaporation::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy,
    YInfType YInf
)
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    const volScalarField& alpha = film.alpha();
    const volScalarField& delta = film.delta();
    const volScalarField& rho = film.rho();
    const surfaceScalarField& phi = film.phi();

    // Set local liquidThermo properties
    const liquidProperties& liquidThermo =
        refCast<const heRhoThermopureMixtureliquidProperties>(film.thermo())
       .cellThermoMixture(0).properties();

    const basicSpecieMixture& primarySpecieThermo =
        refCast<const basicSpecieMixture>(film.primaryThermo());

    // Retrieve fields from film model
    const scalarField& pInf = film.pPrimary();
    const scalarField& T = film.thermo().T();
    const scalarField& he = film.thermo().he();
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& V = film.regionMesh().V();
    const scalarField& magSf = film.magSf();
    const scalarField& VbyA = film.VbyA();
    const vectorField dU(film.UPrimary() - film.Us());
    const scalarField limMass
    (
        max(scalar(0), availableMass - deltaMin_*rho*magSf)
    );

    // Molecular weight of vapour [kg/kmol]
    const scalar Wvap = this->Wvap();

    const scalar Wwax = Wwax_.value();
    const scalar Wsolvent = Wsolvent_.value();

    volScalarField::Internal evapRateCoeff
    (
        IOobject
        (
            IOobject::modelName("evapRateCoeff", typeName),
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        film.regionMesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    );

    volScalarField::Internal evapRateInf
    (
        IOobject
        (
            IOobject::modelName("evapRateInf", typeName),
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        film.regionMesh(),
        dimensionedScalar(evapRateCoeff.dimensions(), 0)
    );

    // Local surface temperature at which evaporation takes place
    scalarField Tloc(dMass.size());

    bool filmPresent = false;

    forAll(dMass, celli)
    {
        if (delta[celli] > deltaMin_)
        {
            filmPresent = true;

            const scalar Ysolvent = Ysolvent_[celli];

            // Molefraction of solvent in liquid film
            const scalar Xsolvent
            (
                Ysolvent*Wsolvent/((1 - Ysolvent)*Wwax + Ysolvent*Wsolvent)
            );

            // Primary region density [kg/m^3]
            const scalar rhoInfc = rhoInf[celli];

            // Cell pressure [Pa]
            const scalar pc = pInf[celli];

            // Calculate the boiling temperature
            const scalar Tb = liquidThermo.pvInvert(pc);

            // Local temperature - impose lower limit of 200 K for stability
            Tloc[celli] = min(TbFactor_*Tb, max(200.0, T[celli]));

            const scalar pPartialCoeff
            (
                liquidThermo.pv(pc, Tloc[celli])*activityCoeff_->value(Xsolvent)
            );

            scalar XsCoeff = pPartialCoeff/pc;

            // Vapour phase mole fraction of solvent at interface
            scalar Xs = XsCoeff*Xsolvent;

            if (Xs > 1)
            {
                WarningInFunction
                    << "Solvent vapour pressure > ambient pressure"
                    << endl;

                XsCoeff /= Xs;
                Xs = 1;
            }

            // Vapour phase mass fraction of solvent at the interface
            const scalar YsCoeff
            (
                XsCoeff/(XsCoeff*Xsolvent*Wsolvent + (1 - Xs)*Wvap)
            );

            // Primary region viscosity [Pa.s]
            const scalar muInfc = muInf[celli];

            // Reynolds number
            const scalar Re = rhoInfc*mag(dU[celli])*L_/muInfc;

            // Vapour diffusivity [m^2/s]
            const scalar Dab = liquidThermo.D(pc, Tloc[celli]);

            // Schmidt number
            const scalar Sc = muInfc/(rhoInfc*(Dab + rootVSmall));

            // Sherwood number
            const scalar Sh = this->Sh(Re, Sc);

            // Mass transfer coefficient [kg/m^3 s]
            evapRateCoeff[celli] =
                rhoInfc*Sh*Dab/(VbyA[celli]*(L_ + rootVSmall));

            // Solvent mass transfer
            const scalar dm
            (
                max
                (
                    dt*V[celli]
                   *evapRateCoeff[celli]*(YsCoeff*Ysolvent - YInf[celli]),
                    0
                )
            );

            if (dm > limMass[celli])
            {
                evapRateCoeff[celli] *= limMass[celli]/dm;
            }

            evapRateInf[celli] = evapRateCoeff[celli]*YInf[celli];
            evapRateCoeff[celli] *= YsCoeff;
        }
    }

    const dimensionedScalar rho0Bydt
    (
        "rho0Bydt",
        dimDensity/dimTime,
        rootVSmall/dt
    );

    volScalarField::Internal impingementRate
    (
        max
        (
           -film.rhoSp(),
            dimensionedScalar(film.rhoSp().dimensions(), 0)
        )
    );

    if (filmPresent)
    {
        // Solve for the solvent mass fraction
        fvScalarMatrix YsolventEqn
        (
            fvm::ddt(alpha, rho, Ysolvent_) + fvm::div(phi, Ysolvent_)
          - fvm::Sp(film.continuityErr(), Ysolvent_)
         ==
            rho0Bydt*Ysolvent_()

          + evapRateInf

            // Include the effect of the impinging droplets
            // added with Ysolvent = Ysolvent0
          + impingementRate*Ysolvent0_

          - fvm::Sp
            (
                rho0Bydt
              + evapRateCoeff
              + film.rhoSp()
              + impingementRate,
                Ysolvent_
            )
        );

        YsolventEqn.relax();
        YsolventEqn.solve();

        Ysolvent_.min(1);
        Ysolvent_.max(0);

        scalarField dm
        (
            dt*V*rhoInf*(evapRateCoeff*Ysolvent_ + evapRateInf)
        );

        dMass += dm;

        // Assume that the vapour transferred to the primary region is
        // already at temperature Tloc so that all heat required for
        // the phase-change is provided by the film
        dEnergy += dm*primarySpecieThermo.Hs(vapId(), pInf, Tloc);

        // Heat is assumed to be removed by heat-transfer to the wall
        // so the energy remains unchanged by the phase-change.
        dEnergy += dm*he;
        // dEnergy += dm*(h[celli] + hVap);
    }
}


void waxSolventEvaporation::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    if (YInfZero_)
    {
        correctModel(dt, availableMass, dMass, dEnergy, zeroField());
    }
    else
    {
        const thermoSingleLayer& film = filmType<thermoSingleLayer>();
        const scalarField& YInf = film.YPrimary()[vapId()];

        correctModel(dt, availableMass, dMass, dEnergy, YInf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
