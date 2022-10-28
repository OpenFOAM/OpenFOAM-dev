/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "momentumSurfaceFilm.H"

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

#include "mappedWallPolyPatch.H"
#include "distributionMap.H"
#include "filmViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(momentumSurfaceFilm, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::momentumSurfaceFilm::read()
{
    return surfaceFilm::read();
}


void Foam::regionModels::momentumSurfaceFilm::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    rhoSpPrimary_ == Zero;
    USpPrimary_ == Zero;
    pSpPrimary_ == Zero;
}


void Foam::regionModels::momentumSurfaceFilm::
transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    // Update fields from primary region via direct mapped
    // (coupled) boundary conditions
    UPrimary_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    rhoPrimary_.correctBoundaryConditions();
    muPrimary_.correctBoundaryConditions();
}


void Foam::regionModels::momentumSurfaceFilm::
transferPrimaryRegionSourceFields()
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


Foam::tmp<Foam::volScalarField>
Foam::regionModels::momentumSurfaceFilm::pc()
{
    return -fvc::laplacian(sigma(), delta_);
}


Foam::tmp<Foam::volScalarField>
Foam::regionModels::momentumSurfaceFilm::pe()
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
        typedName("pe"),
        p_                             // Pressure (mapped from primary region)
      - tpSp                           // Accumulated particle impingement
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::regionModels::momentumSurfaceFilm::rhog() const
{
    return
        fvc::interpolate
        (
            max(nHat() & -g(), dimensionedScalar(g().dimensions(), 0))*VbyA()
        )*fvc::interpolate(rho());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::regionModels::momentumSurfaceFilm::gGradRho() const
{
    return
        fvc::interpolate
        (
            max(nHat() & -g(), dimensionedScalar(g().dimensions(), 0))*VbyA()
        )*fvc::snGrad(rho());
}


void Foam::regionModels::momentumSurfaceFilm::correctCoverage()
{
    coverage_ == pos(delta_ - deltaSmall_);
}


void Foam::regionModels::momentumSurfaceFilm::updateSubmodels()
{
    DebugInFunction << endl;

    // Update ejection model - mass returned is mass available for ejection
    ejection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryMomentumTrans_);

    const volScalarField::Internal rVDt
    (
        1/(time().deltaT()*regionMesh().V())
    );

    // Update mass source field
    rhoSp_ += rVDt*(cloudMassTrans_() + primaryMassTrans_());
    USp_ += rVDt*(cloudMassTrans_()*U_() + primaryMomentumTrans_());

    momentumTransport_->correct();
}


void Foam::regionModels::momentumSurfaceFilm::predictDelta()
{
    DebugInFunction << endl;

    solve(fvm::ddt(rho(), alpha_) + fvc::div(phi_) == -rhoSp_);

    // Bound film volume fraction
    alpha_.max(0);

    delta_ == alpha_*VbyA();

    // Update continuity error caused by the delta_ bounding
    updateContinuityErr();
}


void Foam::regionModels::momentumSurfaceFilm::updateContinuityErr()
{
    continuityErr_ = (fvc::ddt(alpha_, rho()) + fvc::div(phi_))() + rhoSp_;
}


void Foam::regionModels::momentumSurfaceFilm::continuityCheck()
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


Foam::tmp<Foam::volVectorField::Internal>
Foam::regionModels::momentumSurfaceFilm::Uw() const
{
    tmp<volVectorField::Internal> tUw
    (
        volVectorField::Internal::New
        (
            "Uw",
            regionMesh(),
            dimensionedVector(dimVelocity, Zero)
        )
    );

    volVectorField::Internal& Uw = tUw.ref();

    // Push boundary film velocity values into internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<vector>(Uw, pp.faceCells()) =
            U_.boundaryField()[patchi];
    }

    Uw -= nHat()*(Uw_ & nHat());

    return tUw;
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::regionModels::momentumSurfaceFilm::solveMomentum
(
    const volScalarField& pc,
    const volScalarField& pe
)
{
    DebugInFunction << endl;

    // Evaluate viscosity from user-model
    viscosity_->correct(thermo_->p(), thermo_->T());

    const volScalarField::Internal rVDt
    (
        1/(time().deltaT()*regionMesh().V())
    );

    // Momentum equation
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(alpha_, rho(), U_) + fvm::div(phi_, U_)
      - fvm::Sp(continuityErr_, U_)
     ==
      - USp_
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
                      - fvc::interpolate(rho())*(g() & regionMesh().Sf())
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


void Foam::regionModels::momentumSurfaceFilm::solveAlpha
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
    const surfaceScalarField rhof(fvc::interpolate(rho()));
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
              - rhof*(g() & regionMesh().Sf()),
                0
            )
        )
    );

    surfaceScalarField phid
    (
        "phid",
        // constrainFilmField
        // (
        //     rhof
        //    *constrainPhiHbyA(fvc::flux(HbyA) - alpharAUf*phiu, U_, alpha_),
        //     0
        // )
        rhof*constrainPhiHbyA(fvc::flux(HbyA) - alpharAUf*phiu, U_, alpha_)
    );

    const surfaceScalarField ddrhorAUrhogf
    (
        "alphaCoeff",
        alphaf*rhof*alpharAUf*rhogf
    );

    regionMesh().schemes().setFluxRequired(alpha_.name());

    while (pimple_.correctNonOrthogonal())
    {
        // Film thickness equation
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(rho(), alpha_)
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

Foam::regionModels::momentumSurfaceFilm::momentumSurfaceFilm
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    surfaceFilm(modelType, mesh, g, regionType),
    phaseName_(coeffs_.lookupOrDefault("phase", word::null)),
    pimple_(regionMesh()),

    cumulativeContErr_(0),

    deltaSmall_("deltaSmall", dimLength, small),
    maxCo_(solution().lookupOrDefault<scalar>("maxCo", 0)),

    p_
    (
        IOobject
        (
            "p",
            time().timeName(),
            regionMesh()
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    thermo_(rhoThermo::New(regionMesh())),

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

    primaryMomentumTrans_
    (
        IOobject
        (
            "primaryMomentumTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(dimMass*dimVelocity, Zero),
        zeroGradientFvPatchVectorField::typeName
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
            "U", // Must have same name as U to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(dimVelocity, Zero),
        this->mappedFieldAndInternalPatchTypes<vector>()
    ),

    rhoPrimary_
    (
        IOobject
        (
            // Must have same name as rho to enable mapping
            IOobject::groupName("rho", phaseName_),
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
            // Must have same name as rho to enable mapping
            IOobject::groupName("mu", phaseName_),
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure*dimTime, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    viscosity_(surfaceFilmSubModels::viscosityModel::New(*this, coeffs(), mu_)),

    sigma_(Function1<scalar>::New("sigma", coeffs())),

    availableMass_(regionMesh().nCells(), 0),

    ejection_(*this, coeffs_),

    transfer_(*this, coeffs_),

    momentumTransport_
    (
        surfaceFilmSubModels::momentumTransportModel::New(*this, coeffs_)
    ),

    forces_(*this, coeffs_),

    addedMassTotal_(0)
{
    alpha_ == delta_/VbyA();

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctCoverage();

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
            fvc::flux(alpha_*rho()*U_)
        );

        phi_ == phi;
        phiU_ = fvc::flux(U_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::momentumSurfaceFilm::~momentumSurfaceFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::regionModels::momentumSurfaceFilm::sigma() const
{
    tmp<volScalarField> tsigma
    (
        volScalarField::New
        (
            typedName("sigma"),
            regionMesh(),
            dimensionedScalar(dimMass/sqr(dimTime), 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    tsigma.ref().primitiveFieldRef() = sigma_->value(thermo_->T());

    tsigma.ref().correctBoundaryConditions();

    return tsigma;
}


void Foam::regionModels::momentumSurfaceFilm::addSources
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


void Foam::regionModels::momentumSurfaceFilm::preEvolveRegion()
{
    DebugInFunction << endl;

    surfaceFilm::preEvolveRegion();

    transferPrimaryRegionThermoFields();

    transferPrimaryRegionSourceFields();

    // Reset transfer fields
    availableMass_ = mass();
    cloudMassTrans_ == dimensionedScalar(dimMass, 0);
    cloudDiameterTrans_ == dimensionedScalar(dimLength, 0);
    primaryMassTrans_ == dimensionedScalar(dimMass, 0);
    primaryMomentumTrans_ == dimensionedVector(dimMass*dimVelocity, Zero);
}


void Foam::regionModels::momentumSurfaceFilm::evolveRegion()
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


Foam::scalar
Foam::regionModels::momentumSurfaceFilm::CourantNumber() const
{
    const scalarField sumPhi(fvc::surfaceSum(mag(phiU_))().primitiveField());

    return 0.5*gMax(sumPhi/regionMesh().V().field())*time_.deltaTValue();
}


Foam::scalar
Foam::regionModels::momentumSurfaceFilm::maxDeltaT() const
{
    if (maxCo_ > 0)
    {
        return maxCo_*time_.deltaTValue()/(CourantNumber() + small);
    }
    else
    {
        return great;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::regionModels::momentumSurfaceFilm::primaryMassTrans() const
{
    return primaryMassTrans_;
}


const Foam::volScalarField&
Foam::regionModels::momentumSurfaceFilm::cloudMassTrans() const
{
    return cloudMassTrans_;
}


const Foam::volScalarField&
Foam::regionModels::momentumSurfaceFilm::cloudDiameterTrans() const
{
    return cloudDiameterTrans_;
}


Foam::tmp<Foam::volVectorField>
Foam::regionModels::momentumSurfaceFilm::primaryMomentumTrans() const
{
    return primaryMomentumTrans_;
}


void Foam::regionModels::momentumSurfaceFilm::info()
{
    Info<< "\nSurface film: " << type() << endl;

    const scalarField& deltaInternal = delta_;
    const vectorField& Uinternal = U_;
    scalar addedMassTotal = 0;
    outputProperties().readIfPresent("addedMassTotal", addedMassTotal);
    addedMassTotal += returnReduce(addedMassTotal_, sumOp<scalar>());

    Info<< indent << "added mass         = " << addedMassTotal << nl
        << indent << "current mass       = "
        << gSum((delta_*rho()*magSf())()) << nl
        << indent << "min/max(mag(U))    = " << gMin(mag(Uinternal)) << ", "
        << gMax(mag(Uinternal)) << nl
        << indent << "max Courant number = " << CourantNumber() << nl
        << indent << "min/max(delta)     = " << gMin(deltaInternal) << ", "
        << gMax(deltaInternal) << nl
        << indent << "coverage           = "
        << gSum(coverage_.primitiveField()*magSf())/gSum(magSf()) <<  nl;

    ejection_.info(Info);
    transfer_.info(Info);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::regionModels::momentumSurfaceFilm::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        volScalarField::Internal::New
        (
            typedName("Srho"),
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


Foam::tmp<Foam::volScalarField::Internal>
Foam::regionModels::momentumSurfaceFilm::SYi
(
    const label i
) const
{
    return volScalarField::Internal::New
    (
        typedName("SY(" + Foam::name(i) + ")"),
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    );
}


Foam::tmp<Foam::volVectorField::Internal>
Foam::regionModels::momentumSurfaceFilm::SU() const
{
    tmp<volVectorField::Internal> tSU
    (
        volVectorField::Internal::New
        (
            typedName("SU"),
            primaryMesh(),
            dimensionedVector(dimDensity*dimVelocity/dimTime, Zero)
        )
    );

    vectorField& SU = tSU.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        vectorField patchMomentum =
            primaryMomentumTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchMomentum);

        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchIDs()[i]].faceCells();

        forAll(patchMomentum, j)
        {
            SU[cells[j]] += patchMomentum[j]/(V[cells[j]]*dt);
        }
    }

    return tSU;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::regionModels::momentumSurfaceFilm::Sh() const
{
    return volScalarField::Internal::New
    (
        typedName("Sh"),
        primaryMesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


// ************************************************************************* //
