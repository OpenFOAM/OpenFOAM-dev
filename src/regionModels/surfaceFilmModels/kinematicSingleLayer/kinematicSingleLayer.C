/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "kinematicSingleLayer.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvcFlux.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "constrainHbyA.H"

#include "addToRunTimeSelectionTable.H"
#include "mappedWallPolyPatch.H"
#include "mapDistribute.H"
#include "filmThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, kinematicSingleLayer, mesh);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool kinematicSingleLayer::read()
{
    return surfaceFilmRegionModel::read();
}


void kinematicSingleLayer::correctThermoFields()
{
    rho_ == thermo_->rho();
    mu_ == thermo_->mu();
    sigma_ == thermo_->sigma();
}


void kinematicSingleLayer::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    rhoSpPrimary_ == Zero;
    USpPrimary_ == Zero;
    pSpPrimary_ == Zero;
}


void kinematicSingleLayer::transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    // Update fields from primary region via direct mapped
    // (coupled) boundary conditions
    UPrimary_.correctBoundaryConditions();
    pPrimary_.correctBoundaryConditions();
    rhoPrimary_.correctBoundaryConditions();
    muPrimary_.correctBoundaryConditions();
}


void kinematicSingleLayer::transferPrimaryRegionSourceFields()
{
    DebugInFunction << endl;

    volScalarField::Boundary& rhoSpPrimaryBf =
        rhoSpPrimary_.boundaryFieldRef();

    volVectorField::Boundary& USpPrimaryBf =
        USpPrimary_.boundaryFieldRef();

    volScalarField::Boundary& pSpPrimaryBf =
        pSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(rhoSpPrimary_.boundaryField(), patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1/deltaT)
           /primaryMesh().magSf().boundaryField()[patchi]
        );

        rhoSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
        USpPrimaryBf[patchi] *= rpriMagSfdeltaT;
        pSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region
    toRegion(rhoSp_, rhoSpPrimaryBf);
    rhoSp_.field() /= VbyA();
    toRegion(USp_, USpPrimaryBf);
    USp_.field() /= VbyA();
    toRegion(pSp_, pSpPrimaryBf);

    // update addedMassTotal counter
    if (time().writeTime())
    {
        scalar addedMassTotal = 0;
        outputProperties().readIfPresent("addedMassTotal", addedMassTotal);
        addedMassTotal += returnReduce(addedMassTotal_, sumOp<scalar>());
        outputProperties().add("addedMassTotal", addedMassTotal, true);
        addedMassTotal_ = 0;
    }
}


tmp<volScalarField> kinematicSingleLayer::pc()
{
    return -fvc::laplacian(sigma_, delta_);
}


tmp<volScalarField> kinematicSingleLayer::pe()
{
    tmp<volScalarField> tpSp
    (
        volScalarField::New
        (
            "pSp",
            regionMesh(),
            dimensionedScalar(pSp_.dimensions(), 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tpSp.ref().primitiveFieldRef() = pSp_;
    tpSp.ref().correctBoundaryConditions();

    return volScalarField::New
    (
        IOobject::modelName("pe", typeName),
        pPrimary_                      // Pressure (mapped from primary region)
      - tpSp                           // Accumulated particle impingement
    );
}


tmp<surfaceScalarField> kinematicSingleLayer::rhog() const
{
    return
        fvc::interpolate
        (
            max(nHat() & -g_, dimensionedScalar(g_.dimensions(), 0))*VbyA()
        )*fvc::interpolate(rho_);
}


tmp<surfaceScalarField> kinematicSingleLayer::gGradRho() const
{
    return
        fvc::interpolate
        (
            max(nHat() & -g_, dimensionedScalar(g_.dimensions(), 0))*VbyA()
        )*fvc::snGrad(rho_);
}


void kinematicSingleLayer::correctCoverage()
{
    coverage_ == pos(delta_ - deltaSmall_);
}


void kinematicSingleLayer::updateSubmodels()
{
    DebugInFunction << endl;

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, cloudMassTrans_);

    // Update mass source field
    rhoSp_ += cloudMassTrans_/regionMesh().V()/time().deltaT();

    momentumTransport_->correct();
}


void kinematicSingleLayer::predictDelta()
{
    DebugInFunction << endl;

    solve(fvm::ddt(rho_, alpha_) + fvc::div(phi_) == -rhoSp_);

    // Bound film volume fraction
    alpha_.max(0);

    delta_ == alpha_*VbyA();

    // Update continuity error caused by the delta_ bounding
    updateContinuityErr();
}


void kinematicSingleLayer::updateContinuityErr()
{
    continuityErr_ = (fvc::ddt(alpha_, rho_) + fvc::div(phi_))() + rhoSp_;
}


void kinematicSingleLayer::continuityCheck()
{
    const dimensionedScalar totalMass = fvc::domainIntegrate(mass());

    if (totalMass.value() > small)
    {
        const volScalarField::Internal massErr
        (
            time_.deltaT()*magSf()*continuityErr_
        );

        const scalar sumLocalContErr =
        (
            fvc::domainIntegrate(mag(massErr))/totalMass
        ).value();

        const scalar globalContErr =
        (
            fvc::domainIntegrate(massErr)/totalMass
        ).value();

        cumulativeContErr_ += globalContErr;

        Info<< "time step continuity errors: sum local = "
            << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_
            << endl;
    }
}


void kinematicSingleLayer::updateSurfaceVelocities()
{
    // Push boundary film velocity values into internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<vector>(Uw_, pp.faceCells()) =
            U_.boundaryField()[patchi];
    }

    Uw_ -= nHat()*(Uw_ & nHat());

    Us_ = momentumTransport_->Us();
}


tmp<Foam::fvVectorMatrix> kinematicSingleLayer::solveMomentum
(
    const volScalarField& pc,
    const volScalarField& pe
)
{
    DebugInFunction << endl;

    const volScalarField::Internal rVDt
    (
        1/(time().deltaT()*regionMesh().V())
    );

    // Momentum
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(alpha_, rho_, U_) + fvm::div(phi_, U_)
      - fvm::Sp(continuityErr_, U_)
     ==
      - USp_

        // Temporary treatment for the loss of momentum due to mass loss
        // These transfers are not currently included in USp_
      - fvm::Sp(rVDt*(cloudMassTrans_() + primaryMassTrans_()), U_)

      + forces_.correct(U_)
      + momentumTransport_->Su(U_)
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    if (pimple_.momentumPredictor())
    {
        const surfaceScalarField alphaf(fvc::interpolate(alpha_));

        solve
        (
            UEqn
         ==
           -fvc::reconstruct
            (
                constrainFilmField
                (
                    alphaf
                   *(
                        (
                            fvc::snGrad(pe + pc, "snGrad(p)")
                          + gGradRho()*alphaf
                          + rhog()*fvc::snGrad(alpha_)
                        )*regionMesh().magSf()
                      - fvc::interpolate(rho_)*(g_ & regionMesh().Sf())
                    ), 0
                )
            )
        );

        // Remove any patch-normal components of velocity
        U_ -= nHat()*(nHat() & U_);

        U_.correctBoundaryConditions();
    }

    return tUEqn;
}


void kinematicSingleLayer::solveAlpha
(
    const fvVectorMatrix& UEqn,
    const volScalarField& pc,
    const volScalarField& pe
)
{
    DebugInFunction << endl;

    const volScalarField rAU(1/UEqn.A());
    const volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, alpha_));

    const surfaceScalarField alphaf(fvc::interpolate(alpha_));
    const surfaceScalarField rhof(fvc::interpolate(rho_));
    const surfaceScalarField alpharAUf(fvc::interpolate(alpha_*rAU));
    const surfaceScalarField rhogf(rhog());

    const surfaceScalarField phiu
    (
        "phiu",
        (
            constrainFilmField
            (
                (
                    fvc::snGrad(pe + pc, "snGrad(p)")
                  + gGradRho()*alphaf
                )*regionMesh().magSf()
              - rhof*(g_ & regionMesh().Sf()),
                0
            )
        )
    );

    surfaceScalarField phid
    (
        "phid",
        constrainFilmField(rhof*(fvc::flux(HbyA) - alpharAUf*phiu), 0)
    );

    const surfaceScalarField ddrhorAUrhogf
    (
        "alphaCoeff",
        alphaf*rhof*alpharAUf*rhogf
    );

    regionMesh().setFluxRequired(alpha_.name());

    while (pimple_.correctNonOrthogonal())
    {
        // Film thickness equation
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(rho_, alpha_)
          + fvm::div(phid, alpha_)
          - fvm::laplacian(ddrhorAUrhogf, alpha_)
         ==
           -rhoSp_
        );

        alphaEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ == alphaEqn.flux();

            const surfaceScalarField phiGradAlpha
            (
                constrainFilmField
                (
                    rhogf*fvc::snGrad(alpha_)*regionMesh().magSf(),
                    0
                )
            );

            phiU_ = constrainFilmField(phid/rhof - alpharAUf*phiGradAlpha, 0);

            // Update U field
            U_ = HbyA - rAU*fvc::reconstruct(alphaf*(phiu + phiGradAlpha));

            // Remove any patch-normal components of velocity
            U_ -= nHat()*(nHat() & U_);

            U_.correctBoundaryConditions();
        }
    }

    // Bound film volume fraction
    alpha_.max(0);

    delta_ == alpha_*VbyA();

    updateContinuityErr();

    // Continuity check
    continuityCheck();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicSingleLayer::kinematicSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    surfaceFilmRegionModel(modelType, mesh, g, regionType),
    pimple_(regionMesh()),

    cumulativeContErr_(0),

    deltaSmall_("deltaSmall", dimLength, small),
    deltaCoLimit_(solution().lookupOrDefault("deltaCoLimit", 1e-4)),

    rho_
    (
        IOobject
        (
            "rho",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    mu_
    (
        IOobject
        (
            "mu",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure*dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    sigma_
    (
        IOobject
        (
            "sigma",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    delta_
    (
        IOobject
        (
            "delta",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    alpha_
    (
        IOobject
        (
            "alpha",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        delta_/VbyA(),
        delta_.boundaryField().types()
    ),

    U_
    (
        IOobject
        (
            "U",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    Us_
    (
        IOobject
        (
            "Us",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
    ),

    Uw_
    (
        IOobject
        (
            "Uw",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
    ),

    phi_
    (
        IOobject
        (
            "phi",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime, 0)
    ),

    phiU_
    (
        IOobject
        (
            "phiU",
            time().timeName(),
            regionMesh()
        ),
        regionMesh(),
        dimensionedScalar(dimVolume/dimTime, 0)
    ),

    continuityErr_
    (
        IOobject
        (
            "continuityErr",
            time().timeName(),
            regionMesh()
        ),
        regionMesh(),
        dimensionedScalar(alpha_.dimensions()*dimDensity/dimTime, 0)
    ),

    coverage_
    (
        IOobject
        (
            "coverage",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    primaryMassTrans_
    (
        IOobject
        (
            "primaryMassTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    cloudMassTrans_
    (
        IOobject
        (
            "cloudMassTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    cloudDiameterTrans_
    (
        IOobject
        (
            "cloudDiameterTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimLength, -1),
        zeroGradientFvPatchScalarField::typeName
    ),

    rhoSp_
    (
        IOobject
        (
            "rhoSp",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    ),

    USp_
    (
        IOobject
        (
            "USp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(dimDensity*dimVelocity/dimTime, Zero)
    ),

    pSp_
    (
        IOobject
        (
            "pSp",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, 0)
    ),

    rhoSpPrimary_
    (
        IOobject
        (
            rhoSp_.name(),
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimMass, 0)
    ),

    USpPrimary_
    (
        IOobject
        (
            USp_.name(),
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedVector(dimMass*dimVelocity, Zero)
    ),

    pSpPrimary_
    (
        IOobject
        (
            pSp_.name(),
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimMass*dimVelocity, 0)
    ),

    UPrimary_
    (
        IOobject
        (
            "U", // must have same name as U to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(dimVelocity, Zero),
        this->mappedFieldAndInternalPatchTypes<vector>()
    ),

    pPrimary_
    (
        IOobject
        (
            "p", // must have same name as p to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    rhoPrimary_
    (
        IOobject
        (
            "thermo:rho", // must have same name as rho to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    muPrimary_
    (
        IOobject
        (
            "thermo:mu", // must have same name as mu to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure*dimTime, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    thermo_(thermoModel::New(*this, coeffs_)),

    availableMass_(regionMesh().nCells(), 0),

    injection_(*this, coeffs_),

    transfer_(*this, coeffs_),

    momentumTransport_(momentumTransportModel::New(*this, coeffs_)),

    forces_(*this, coeffs_),

    addedMassTotal_(0)
{
    alpha_ == delta_/VbyA();

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctCoverage();

        correctThermoFields();

        surfaceScalarField phi
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
            fvc::flux(alpha_*rho_*U_)
        );

        phi_ == phi;
        phiU_ = fvc::flux(U_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kinematicSingleLayer::~kinematicSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    DebugInFunction
        << "\nSurface film: " << type() << ": adding to film source:" << nl
        << "    mass     = " << massSource << nl
        << "    momentum = " << momentumSource << nl
        << "    pressure = " << pressureSource << endl;

    rhoSpPrimary_.boundaryFieldRef()[patchi][facei] -= massSource;
    USpPrimary_.boundaryFieldRef()[patchi][facei] -= momentumSource;
    pSpPrimary_.boundaryFieldRef()[patchi][facei] -= pressureSource;

    addedMassTotal_ += massSource;
}


void kinematicSingleLayer::preEvolveRegion()
{
    DebugInFunction << endl;

    surfaceFilmRegionModel::preEvolveRegion();

    transferPrimaryRegionThermoFields();

    correctThermoFields();

    transferPrimaryRegionSourceFields();

    // Reset transfer fields
    availableMass_ = mass();
    cloudMassTrans_ == dimensionedScalar(dimMass, 0);
    cloudDiameterTrans_ == dimensionedScalar(dimLength, 0);
    primaryMassTrans_ == dimensionedScalar(dimMass, 0);
}


void kinematicSingleLayer::evolveRegion()
{
    DebugInFunction << endl;

    // Update film coverage indicator
    correctCoverage();

    // Update film wall and surface velocities
    updateSurfaceVelocities();

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

        // Solve for momentum
        const fvVectorMatrix UEqn(solveMomentum(pc, pe));

        // Film thickness correction loop
        while (pimple_.correct())
        {
            solveAlpha(UEqn, pc, pe);
        }
    }

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


scalar kinematicSingleLayer::CourantNumber() const
{
    const scalarField sumPhi(fvc::surfaceSum(mag(phiU_))().primitiveField());

    const scalar CoNum =
        0.5*gMax(sumPhi/regionMesh().V().field())*time_.deltaTValue();

    Info<< "Film max Courant number: " << CoNum << endl;

    return CoNum;
}


tmp<volScalarField> kinematicSingleLayer::primaryMassTrans() const
{
    return primaryMassTrans_;
}


const volScalarField& kinematicSingleLayer::cloudMassTrans() const
{
    return cloudMassTrans_;
}


const volScalarField& kinematicSingleLayer::cloudDiameterTrans() const
{
    return cloudDiameterTrans_;
}


void kinematicSingleLayer::info()
{
    Info<< "\nSurface film: " << type() << endl;

    const scalarField& deltaInternal = delta_;
    const vectorField& Uinternal = U_;
    scalar addedMassTotal = 0;
    outputProperties().readIfPresent("addedMassTotal", addedMassTotal);
    addedMassTotal += returnReduce(addedMassTotal_, sumOp<scalar>());

    Info<< indent << "added mass         = " << addedMassTotal << nl
        << indent << "current mass       = "
        << gSum((delta_*rho_*magSf())()) << nl
        << indent << "min/max(mag(U))    = " << gMin(mag(Uinternal)) << ", "
        << gMax(mag(Uinternal)) << nl
        << indent << "min/max(delta)     = " << gMin(deltaInternal) << ", "
        << gMax(deltaInternal) << nl
        << indent << "coverage           = "
        << gSum(coverage_.primitiveField()*magSf())/gSum(magSf()) <<  nl;

    injection_.info(Info);
    transfer_.info(Info);
}


tmp<volScalarField::Internal> kinematicSingleLayer::Srho() const
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Srho", typeName),
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    );
}


tmp<volScalarField::Internal> kinematicSingleLayer::Srho
(
    const label i
) const
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Srho(" + Foam::name(i) + ")", typeName),
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    );
}


tmp<volScalarField::Internal> kinematicSingleLayer::Sh() const
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Sh", typeName),
        primaryMesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
