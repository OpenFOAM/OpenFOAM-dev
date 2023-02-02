/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "thermoSurfaceFilm.H"
#include "basicSpecieMixture.H"
#include "liquidThermo.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoSurfaceFilm, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::thermoSurfaceFilm::read()
{
    // No additional properties to read
    return momentumSurfaceFilm::read();
}


void Foam::thermoSurfaceFilm::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    momentumSurfaceFilm::resetPrimaryRegionSourceTerms();

    hSpPrimary_ == dimensionedScalar(hSp_.dimensions(), 0);
}


void Foam::thermoSurfaceFilm::transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    momentumSurfaceFilm::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void Foam::thermoSurfaceFilm::transferPrimaryRegionSourceFields()
{
    DebugInFunction << endl;

    momentumSurfaceFilm::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hSpPrimaryBf = hSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region
    toRegion(hSp_, hSpPrimaryBf);
    hSp_.field() /= VbyA();
}


void Foam::thermoSurfaceFilm::correctCoverage()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(coverage_, i)
        {
            if ((coverage_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                coverage_[i] = 1;
            }
            else if ((coverage_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                coverage_[i] = 0;
            }
        }

        coverage_.correctBoundaryConditions();
    }
    else
    {
        coverage_ == pos(delta_ - dimensionedScalar(dimLength, deltaWet_));
    }
}


void Foam::thermoSurfaceFilm::updateSubmodels()
{
    DebugInFunction << endl;

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update radiation
    radiation_->correct();

    // Update ejection model - mass returned is mass available for ejection
    ejection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    const volScalarField::Internal rMagSfDt((1/time().deltaT())/magSf());

    // Vapour recoil pressure
    pSp_ -= sqr(rMagSfDt*primaryMassTrans_())/(2*rhoPrimary_());

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct
    (
        availableMass_,
        primaryMassTrans_,
        primaryMomentumTrans_,
        primaryEnergyTrans_
    );

    const volScalarField::Internal rVDt
    (
        1/(time().deltaT()*mesh().V())
    );

    volScalarField& he = thermo_->he();

    // Update source fields
    rhoSp_ += rVDt*(cloudMassTrans_() + primaryMassTrans_());
    USp_ += rVDt*(cloudMassTrans_()*U_() + primaryMomentumTrans_());
    hSp_ += rVDt*(cloudMassTrans_()*he() + primaryEnergyTrans_());

    momentumTransport_->correct();
}


Foam::tmp<Foam::fvScalarMatrix> Foam::thermoSurfaceFilm::q
(
    volScalarField& h
) const
{
    const volScalarField::Internal coverage(pos(delta_() - deltaSmall_));

    const volScalarField::Internal& T = thermo().T();

    const tmp<volScalarField> tCpv = thermo().Cpv();
    const volScalarField::Internal& Cpv = tCpv();

    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp((htcs_->h()/VbyA())/Cpv, h)
      + (htcs_->h()/VbyA())*(h()/Cpv + coverage*(TPrimary_() - T))

        // Heat-transfer to the wall
      - fvm::Sp((htcw_->h()/VbyA())/Cpv, h)
      + (htcw_->h()/VbyA())*(h()/Cpv + coverage*(Tw() - T))
    );
}


void Foam::thermoSurfaceFilm::solveEnergy()
{
    DebugInFunction << endl;

    volScalarField& he = thermo_->he();

    fvScalarMatrix heEqn
    (
        fvm::ddt(alpha_, rho(), he) + fvm::div(phi_, he)
      - fvm::Sp(continuityErr_, he)
     ==
      - hSp_
      + q(he)
      + radiation_->Shs()/VbyA()
    );

    heEqn.relax();

    heEqn.solve();

    thermo_->correct();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoSurfaceFilm::thermoSurfaceFilm
(
    const word& modelType,
    const fvMesh& primaryMesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    momentumSurfaceFilm(modelType, primaryMesh, g, regionType, false),

    primaryThermo_
    (
        primaryMesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        )
    ),

    primaryEnergyTrans_
    (
        IOobject
        (
            "primaryEnergyTrans",
            time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    deltaWet_(coeffs_.lookup<scalar>("deltaWet")),
    hydrophilic_(readBool(coeffs_.lookup("hydrophilic"))),
    hydrophilicDryScale_(0),
    hydrophilicWetScale_(0),

    hSp_
    (
        IOobject
        (
            "hSp",
            time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    ),

    hSpPrimary_
    (
        IOobject
        (
            hSp_.name(),
            time().name(),
            primaryMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh,
        dimensionedScalar(hSp_.dimensions(), 0)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // Same name as T on primary region to enable mapping
            time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimTemperature, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    YPrimary_(),

    htcs_
    (
        surfaceFilmModels::heatTransferModel::New
        (
            *this,
            coeffs().subDict("upperSurfaceModels")
        )
    ),

    htcw_
    (
        surfaceFilmModels::heatTransferModel::New
        (
            *this,
            coeffs().subDict("lowerSurfaceModels")
        )
    ),

    phaseChange_(surfaceFilmModels::phaseChangeModel::New(*this, coeffs())),
    radiation_(surfaceFilmModels::radiationModel::New(*this, coeffs())),
    Tmin_(-vGreat),
    Tmax_(vGreat)
{
    if (coeffs().readIfPresent("Tmin", Tmin_))
    {
        Info<< "    limiting minimum temperature to " << Tmin_ << endl;
    }

    if (coeffs().readIfPresent("Tmax", Tmax_))
    {
        Info<< "    limiting maximum temperature to " << Tmax_ << endl;
    }

    if (isA<basicSpecieMixture>(primaryThermo_))
    {
        const basicSpecieMixture& primarySpecieThermo =
            refCast<const basicSpecieMixture>(primaryThermo_);

        YPrimary_.setSize(primarySpecieThermo.species().size());

        forAll(primarySpecieThermo.species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        primarySpecieThermo.species()[i],
                        time().name(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar(dimless, 0),
                    this->mappedFieldAndInternalPatchTypes<scalar>()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        coeffs_.lookup("hydrophilicDryScale") >> hydrophilicDryScale_;
        coeffs_.lookup("hydrophilicWetScale") >> hydrophilicWetScale_;
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctCoverage();

        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                time().name(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            fvc::flux(alpha_*rho()*U_)
        );

        phi_ == phi;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoSurfaceFilm::~thermoSurfaceFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermoSurfaceFilm::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    momentumSurfaceFilm::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    DebugInFunction << "    energy   = " << energySource << endl;

    hSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void Foam::thermoSurfaceFilm::preEvolveRegion()
{
    DebugInFunction << endl;

    momentumSurfaceFilm::preEvolveRegion();
    primaryEnergyTrans_ == dimensionedScalar(dimEnergy, 0);
}


void Foam::thermoSurfaceFilm::evolveRegion()
{
    DebugInFunction << endl;

    // Update film coverage indicator
    correctCoverage();

    // Predict delta_ from continuity
    predictDelta();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Predict delta_ from continuity with updated source
    predictDelta();

    // Capillary pressure
    const volScalarField pc(this->pc());

    while (pimple_.loop())
    {
        // External pressure
        const volScalarField pe(this->pe());

        // Solve for momentum for U_
        const fvVectorMatrix UEqn(solveMomentum(pc, pe));

        // Solve energy for h_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        while (pimple_.correct())
        {
            solveAlpha(UEqn, pc, pe);
        }
    }

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::thermoSurfaceFilm::Ts() const
{
    return thermo().T();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::thermoSurfaceFilm::Tw() const
{
    tmp<volScalarField::Internal> tTw
    (
        volScalarField::Internal::New
        (
            "Tw",
            mesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    volScalarField::Internal& Tw = tTw.ref();

    const volScalarField& T = thermo().T();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        UIndirectList<scalar>(Tw, pp.faceCells()) =
            T.boundaryField()[patchi];
    }

    return tTw;
}


void Foam::thermoSurfaceFilm::info()
{
    momentumSurfaceFilm::info();

    const scalarField& Tinternal = thermo().T();

    Info<< indent << "min/mean/max(T)    = "
        << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", "
        << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::thermoSurfaceFilm::SYi(const label i) const
{
    const basicSpecieMixture& primarySpecieThermo =
        refCast<const basicSpecieMixture>(primaryThermo_);

    // Set local liquidThermo properties
    const liquidProperties& liquidThermo =
        refCast<const heRhoThermopureMixtureliquidProperties>(thermo())
       .cellThermoMixture(0).properties();

    const label vapId = primarySpecieThermo.species()[liquidThermo.name()];

    tmp<volScalarField::Internal> tSYi
    (
        volScalarField::Internal::New
        (
            typedName("SY(" + Foam::name(i) + ")"),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    if (vapId == i)
    {
        scalarField& SYi = tSYi.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                SYi[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSYi;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::thermoSurfaceFilm::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        volScalarField::Internal::New
        (
            typedName("Sh"),
            primaryMesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    scalarField& Sh = tSh.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchIDs()[i]].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }

    return tSh;
}


// ************************************************************************* //
