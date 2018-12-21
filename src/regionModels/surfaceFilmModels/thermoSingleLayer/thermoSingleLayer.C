/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "thermoSingleLayer.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcFlux.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, thermoSingleLayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList thermoSingleLayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchi)
    {
        if
        (
            T_.boundaryField()[patchi].fixesValue()
         || isA<mixedFvPatchScalarField>(T_.boundaryField()[patchi])
         || isA<mappedFieldFvPatchField<scalar>>(T_.boundaryField()[patchi])
        )
        {
            bTypes[patchi] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar(hsSp_.dimensions(), 0);
}


void thermoSingleLayer::correctThermoFields()
{
    rho_ == filmThermo_->rho();
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();
}


void thermoSingleLayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    volScalarField::Boundary& hsBf = hs_.boundaryFieldRef();

    forAll(hsBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchi];
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            hsBf[patchi] == hs(Tp, patchi);
        }
    }
}


void thermoSingleLayer::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<scalar>(Tw_, pp.faceCells()) =
            T_.boundaryField()[patchi];
    }
    Tw_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hsSpPrimaryBf =
        hsSpPrimary_.boundaryFieldRef();

    // Convert accummulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hsSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();
}


void thermoSingleLayer::correctAlpha()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos0(delta_ - dimensionedScalar(dimLength, deltaWet_));
    }
}


void thermoSingleLayer::updateSubmodels()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update radiation
    radiation_->correct();

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    const volScalarField rMagSfDt((1/time().deltaT())/magSf());

    // Vapour recoil pressure
    pSp_ -= sqr(rMagSfDt*primaryMassTrans_)/(2*rhoPrimary_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryEnergyTrans_);

    // Update source fields
    rhoSp_ += rMagSfDt*(cloudMassTrans_ + primaryMassTrans_);
    hsSp_ += rMagSfDt*(cloudMassTrans_*hs_ + primaryEnergyTrans_);

    turbulence_->correct();
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& hs) const
{
    const volScalarField alpha(pos(delta_ - deltaSmall_));

    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp(htcs_->h()/Cp_, hs)
      + htcs_->h()*(hs/Cp_ + alpha*(TPrimary_ - T_))

        // Heat-transfer to the wall
      - fvm::Sp(htcw_->h()/Cp_, hs)
      + htcw_->h()*(hs/Cp_ + alpha*(Tw_- T_))
    );
}


void thermoSingleLayer::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    updateSurfaceTemperatures();

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
      - hsSp_
      + q(hs_)
      + radiation_->Shs()
    );

    correctThermoFields();

    // Evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, regionType, false),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    Cp_
    (
        IOobject
        (
            "Cp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimLength/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    T_
    (
        IOobject
        (
            "Tf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    hs_
    (
        IOobject
        (
            "hf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass, 0),
        hsBoundaryTypes()
    ),

    primaryEnergyTrans_
    (
        IOobject
        (
            "primaryEnergyTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    deltaWet_(readScalar(coeffs_.lookup("deltaWet"))),
    hydrophilic_(readBool(coeffs_.lookup("hydrophilic"))),
    hydrophilicDryScale_(0.0),
    hydrophilicWetScale_(0.0),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimArea/dimTime, 0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // Must have same name as hSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(hsSp_.dimensions(), 0)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    YPrimary_(),

    viscosity_(filmViscosityModel::New(*this, coeffs(), mu_)),
    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    radiation_(filmRadiationModel::New(*this, coeffs())),
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

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, 0),
                    pSp_.boundaryField().types()
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

        correctAlpha();

        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_);

        deltaRho_ == delta_*rho_;

        surfaceScalarField phi0
        (
            IOobject
            (
                "phi",
                time().timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            fvc::flux(deltaRho_*U_)
        );

        phi_ == phi0;

        // Evaluate viscosity from user-model
        viscosity_->correct(pPrimary_, T_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        Info<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void thermoSingleLayer::preEvolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::preEvolveRegion();
    primaryEnergyTrans_ == dimensionedScalar(dimEnergy, 0);
}


void thermoSingleLayer::evolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Update film coverage indicator
    correctAlpha();

    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Update film wall and surface temperatures
    updateSurfaceTemperatures();

    // Solve continuity for deltaRho_
    solveContinuity();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update temperature using latest hs_
    T_ == T(hs_);

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const volScalarField& thermoSingleLayer::Cp() const
{
    return Cp_;
}


const volScalarField& thermoSingleLayer::kappa() const
{
    return kappa_;
}


const volScalarField& thermoSingleLayer::T() const
{
    return T_;
}


const volScalarField& thermoSingleLayer::Ts() const
{
    return Ts_;
}


const volScalarField& thermoSingleLayer::Tw() const
{
    return Tw_;
}


const volScalarField& thermoSingleLayer::hs() const
{
    return hs_;
}


void thermoSingleLayer::info()
{
    kinematicSingleLayer::info();

    const scalarField& Tinternal = T_;

    Info<< indent << "min/mean/max(T)    = "
        << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", "
        << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        volScalarField::Internal::New
        (
            "thermoSingleLayer::Srho",
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    scalarField& Srho = tSrho.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchi = intCoupledPatchIDs()[i];

        scalarField patchMass =
            primaryMassTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchMass);

        const label primaryPatchi = primaryPatchIDs()[i];
        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    tmp<volScalarField::Internal> tSrho
    (
        volScalarField::Internal::New
        (
            "thermoSingleLayer::Srho(" + Foam::name(i) + ")",
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const unallocLabelList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        volScalarField::Internal::New
        (
            "thermoSingleLayer::Sh",
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

        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchIDs()[i]].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }

    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
