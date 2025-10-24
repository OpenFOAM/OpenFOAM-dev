/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "phaseSurfaceBoiling.H"

#include "multiphaseEuler.H"

#include "diameterModel.H"
#include "twoResistanceHeatTransfer.H"

#include "saturationTemperatureModel.H"
#include "partitioningModel.H"
#include "nucleationSiteModel.H"
#include "departureDiameterModel.H"
#include "departureFrequencyModel.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseSurfaceBoiling, 0);
    addToRunTimeSelectionTable(fvModel, phaseSurfaceBoiling, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseSurfaceBoiling::readCoeffs(const dictionary& dict)
{
    saturationModelPtr_.reset
    (
        saturationTemperatureModel::New
        (
            "saturationTemperature",
            dict
        ).ptr()
    );

    partitioningModel_.reset
    (
        wallBoilingModels::partitioningModel::New
        (
            dict.subDict("partitioningModel")
        ).ptr()
    );

    nucleationSiteModel_.reset
    (
        wallBoilingModels::nucleationSiteModel::New
        (
            dict.subDict("nucleationSiteModel")
        ).ptr()
    );

    departureDiameterModel_.reset
    (
        wallBoilingModels::departureDiameterModel::New
        (
            dict.subDict("departureDiameterModel")
        ).ptr()
    );

    departureFrequencyModel_.reset
    (
        wallBoilingModels::departureFrequencyModel::New
        (
            dict.subDict("departureFrequencyModel")
        ).ptr()
    );
}


void Foam::fv::phaseSurfaceBoiling::correctMDot() const
{
    using constant::mathematical::pi;

    static const dimensionedScalar rooVSmallT(dimTemperature, rootVSmall);
    static const dimensionedScalar rootVSmallF(dimless/dimTime, rootVSmall);
    static const dimensionedScalar rootVSmallH
    (
        heatTransferModel::dimK,
        rootVSmall
    );

    const rhoThermo& solidThermo = solid_.thermo();
    const volScalarField::Internal& solidT = solidThermo.T();
    const rhoThermo& liquidThermo = liquid_.thermo();
    const volScalarField::Internal& liquidT = liquidThermo.T();
    const rhoThermo& vapourThermo = vapour_.thermo();

    // Estimate the surface temperature from the surrounding temperature and
    // heat transfer coefficients. Note that lagged values of qEvaporative and
    // qQuenching are used in this calculation.
    const Pair<tmp<volScalarField>> Hs =
        solver_.heatTransfer.Hs(liquid_, solid_, scalar(0));
    const volScalarField::Internal& liquidH = Hs.first();
    const volScalarField::Internal& solidH = Hs.second();
    const volScalarField::Internal Tsurface
    (
        (solidH*solidT + liquidH*liquidT + qEvaporative_ + qQuenching_)
       /max(solidH + liquidH, rootVSmallH)
    );

    const volScalarField::Internal Tsat
    (
        saturationModelPtr_->Tsat(liquid_.fluidThermo().p()())
    );

    const volScalarField::Internal L
    (
        vapourThermo.ha(liquid_.fluidThermo().p()(), Tsat)
      - liquidThermo.ha()()
    );

    // Wetted fraction
    wetFraction_ =
        partitioningModel_->wetFraction(liquid_()/max(1 - solid_(), small));

    // Bubble departure diameter
    dDeparture_ =
        departureDiameterModel_->dDeparture
        (
            liquid_,
            vapour_,
            solid_,
            Tsurface,
            Tsat,
            L
        );

    // Bubble departure frequency
    fDeparture_ =
        departureFrequencyModel_->fDeparture
        (
            liquid_,
            vapour_,
            solid_,
            Tsurface,
            Tsat,
            L,
            dDeparture_
        );

    // Nucleation site density
    nucleationSiteDensity_ =
        nucleationSiteModel_->nucleationSiteDensity
        (
            liquid_,
            vapour_,
            solid_,
            Tsurface,
            Tsat,
            L,
            dDeparture_,
            fDeparture_
        );

    const tmp<volScalarField> tliquidRho(liquidThermo.rho());
    const volScalarField::Internal& liquidRho = tliquidRho();
    const tmp<volScalarField> tvapourRho(vapourThermo.rho());
    const volScalarField::Internal& vapourRho = tvapourRho();
    const volScalarField::Internal& liquidCp = liquidThermo.Cp();
    const volScalarField::Internal& liquidkappa = liquidThermo.kappa();

    // Area fractions: Del Valle & Kenning (1985)
    const volScalarField::Internal Ja
    (
        liquidRho*liquidCp*(Tsat - min(Tsurface, Tsat))/(liquidRho*L)
    );
    const volScalarField::Internal Al
    (
        wetFraction_*4.8*exp(min(-Ja/80, log(vGreat)))
    );
    const volScalarField::Internal A2
    (
        min(pi*sqr(dDeparture_)*nucleationSiteDensity_*Al/4, scalar(1))
    );
    const volScalarField::Internal A2E
    (
        min(pi*sqr(dDeparture_)*nucleationSiteDensity_*Al/4, scalar(5))
    );

    // Surface area per unit volume
    const volScalarField Av(solid_.diameter().Av());

    // Relaxation factor
    const scalar f = mesh().solution().fieldRelaxationFactor(mDot_.member());

    // Volumetric mass source in due to the wall boiling and bulk nucleation
    mDot_ =
        (1 - f)*mDot_
      + f*(1.0/6.0)*A2E*dDeparture_*vapourRho*fDeparture_*Av;

    // Evaporative heat flux
    qEvaporative_ = mDot_*L;

    // Quenching heat transfer coefficient
    const volScalarField::Internal hQuenching
    (
        2
       *liquidkappa
       *fDeparture_
       *sqrt
        (
            0.8
           /max(fDeparture_, rootVSmallF)
           /(pi*(liquidkappa/liquidCp)/liquidRho)
        )
    );

    // Quenching heat flux
    qQuenching_ =
        (1 - f)*qQuenching_
      + f*A2*hQuenching*max(Tsurface - liquidThermo.T(), rooVSmallT)*Av;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseSurfaceBoiling::phaseSurfaceBoiling
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange(name, modelType, mesh, dict, wordList()),
    nucleation(),
    solver_(mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName)),
    fluid_(solver_.fluid),
    liquid_(fluid_.phases()[phaseNames().first()]),
    vapour_(fluid_.phases()[phaseNames().second()]),
    solid_(fluid_.phases()[dict.lookup("phase")]),
    saturationModelPtr_(nullptr),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiameterModel_(nullptr),
    departureFrequencyModel_(nullptr),
    pressureEquationIndex_(-1),
    wetFraction_
    (
        IOobject
        (
            name + ":wetFraction",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, scalar(0))
    ),
    dDeparture_
    (
        IOobject
        (
            name + ":dDeparture",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, vGreat)
    ),
    fDeparture_
    (
        IOobject
        (
            name + ":fDeparture",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(inv(dimTime), scalar(0))
    ),
    nucleationSiteDensity_
    (
        IOobject
        (
            name + ":nucleationSiteDensity",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(inv(dimArea), scalar(0))
    ),
    qQuenching_
    (
        IOobject
        (
            name + ":qQuenching",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, scalar(0))
    ),
    qEvaporative_
    (
        IOobject
        (
            name + ":qEvaporative",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, scalar(0))
    ),
    mDot_
    (
        IOobject
        (
            name + ":mDot",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, scalar(0))
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::phaseSurfaceBoiling::addsSupToField(const word& fieldName) const
{
    return
        phaseChange::addsSupToField(fieldName)
     || fieldName == solid_.thermo().he().name();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceBoiling::Lfraction() const
{
    // Put all the latent heat into the liquid
    return
        volScalarField::Internal::New
        (
            name() + ":Lfraction",
            mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceBoiling::d() const
{
    return dDeparture_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceBoiling::nDot() const
{
    return fDeparture_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceBoiling::tau() const
{
    NotImplemented;
    return tmp<volScalarField::Internal>(nullptr);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceBoiling::mDot() const
{
    return mDot_;
}


void Foam::fv::phaseSurfaceBoiling::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &liquid_ || &alpha == &vapour_)
     && (&rho == &liquid_.rho() || &rho == &vapour_.rho())
     && &eqn.psi() == &solver_.p_rgh
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // the current phase loop
        if (pressureEquationIndex_ % 2 == 0) correctMDot();
        pressureEquationIndex_ ++;
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


void Foam::fv::phaseSurfaceBoiling::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    // Liquid or solid energy equation. Apply the additional explicit latent
    // and quenching heat transfers.
    if (&he == &liquid_.thermo().he() || &he == &solid_.thermo().he())
    {
        const scalar sign = &he == &solid_.thermo().he() ? -1 : +1;

        eqn += sign*(qEvaporative_ + qQuenching_);
    }

    // Let the base class do the other liquid-vapour transfers
    if (&he != &solid_.thermo().he())
    {
        phaseChange::addSup(alpha, rho, he, eqn);
    }
}


void Foam::fv::phaseSurfaceBoiling::correct()
{
    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctMDot();
}


bool Foam::fv::phaseSurfaceBoiling::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
