/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "phaseSystem.H"
#include "interfaceSurfaceTensionModel.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"
#include "localEulerDdtScheme.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvCorrectPhi.H"
#include "fvcMeshPhi.H"
#include "correctContactAngle.H"
#include "fixedValueFvsPatchFields.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
    defineRunTimeSelectionTable(phaseSystem, dictionary);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::sumAlphaMoving() const
{
    tmp<volScalarField> sumAlphaMoving
    (
        volScalarField::New
        (
            "sumAlphaMoving",
            movingPhaseModels_[0],
            calculatedFvPatchScalarField::typeName
        )
    );

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        sumAlphaMoving.ref() += movingPhaseModels_[movingPhasei];
    }

    return sumAlphaMoving;
}


void Foam::phaseSystem::setMixtureU(const volVectorField& Um0)
{
    // Calculate the mean velocity difference with respect to Um0
    // from the current velocity of the moving phases
    volVectorField dUm(Um0);

    forAll(movingPhaseModels_, movingPhasei)
    {
        dUm -=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].U();
    }

    forAll(movingPhaseModels_, movingPhasei)
    {
        movingPhaseModels_[movingPhasei].URef() += dUm;
    }
}


void Foam::phaseSystem::setMixturePhi
(
    const PtrList<surfaceScalarField>& alphafs,
    const surfaceScalarField& phim0
)
{
    // Calculate the mean flux difference with respect to phim0
    // from the current flux of the moving phases
    surfaceScalarField dphim(phim0);

    forAll(movingPhaseModels_, movingPhasei)
    {
        dphim -=
            alphafs[movingPhaseModels_[movingPhasei].index()]
           *movingPhaseModels_[movingPhasei].phi();
    }

    forAll(movingPhaseModels_, movingPhasei)
    {
        movingPhaseModels_[movingPhasei].phiRef() += dphim;
    }
}


Foam::tmp<Foam::surfaceVectorField> Foam::phaseSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::K
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(phase1, phase2);

    correctContactAngle
    (
        phase1,
        phase2,
        phase1.U()().boundaryField(),
        deltaN_,
        tnHatfv.ref().boundaryFieldRef()
    );

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    pimple_(mesh_.lookupObject<pimpleNoLoopControl>("solutionControl")),

    MRF_(mesh_),

    referencePhaseName_(lookupOrDefault("referencePhase", word::null)),

    phaseModels_
    (
        lookup("phases"),
        phaseModel::iNew(*this, referencePhaseName_)
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimFlux, 0)
    ),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure/dimTime, 0)
    ),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    // Groupings
    label movingPhasei = 0;
    label stationaryPhasei = 0;
    label anisothermalPhasei = 0;
    label multicomponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        movingPhasei += !phase.stationary();
        stationaryPhasei += phase.stationary();
        anisothermalPhasei += !phase.isothermal();
        multicomponentPhasei += !phase.pure();
    }
    movingPhaseModels_.resize(movingPhasei);
    stationaryPhaseModels_.resize(stationaryPhasei);
    anisothermalPhaseModels_.resize(anisothermalPhasei);
    multicomponentPhaseModels_.resize(multicomponentPhasei);

    movingPhasei = 0;
    stationaryPhasei = 0;
    anisothermalPhasei = 0;
    multicomponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        if (!phase.stationary())
        {
            movingPhaseModels_.set(movingPhasei++, &phase);
        }
        if (phase.stationary())
        {
            stationaryPhaseModels_.set(stationaryPhasei++, &phase);
        }
        if (!phase.isothermal())
        {
            anisothermalPhaseModels_.set(anisothermalPhasei++, &phase);
        }
        if (!phase.pure())
        {
            multicomponentPhaseModels_.set(multicomponentPhasei++, &phase);
        }
    }

    forAll(movingPhaseModels_, movingPhasei)
    {
        phi_ +=
            fvc::interpolate(movingPhaseModels_[movingPhasei])
           *movingPhaseModels_[movingPhasei].phi();
    }

    // Write phi
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    // Interface compression coefficients
    if (this->found("interfaceCompression"))
    {
        generateInterfacialValues("interfaceCompression", cAlphas_);
    }

    // Surface tension models
    generateInterfacialModels(interfaceSurfaceTensionModels_);

    // Update motion fields
    correctKinematics();

    // Set the optional reference phase fraction from the other phases
    if (referencePhaseName_ != word::null)
    {
        phaseModel* referencePhasePtr = &phases()[referencePhaseName_];
        volScalarField& referenceAlpha = *referencePhasePtr;

        referenceAlpha = 1;

        forAll(phaseModels_, phasei)
        {
            if (&phaseModels_[phasei] != referencePhasePtr)
            {
                referenceAlpha -= phaseModels_[phasei];
            }
        }
    }

    forAll(phaseModels_, phasei)
    {
        const volScalarField& alphai = phases()[phasei];
        mesh_.schemes().setFluxRequired(alphai.name());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    tmp<volScalarField> rho(movingPhaseModels_[0]*movingPhaseModels_[0].rho());

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        rho.ref() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].rho();
    }

    if (stationaryPhaseModels_.empty())
    {
        return rho;
    }
    else
    {
        return rho/sumAlphaMoving();
    }
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    tmp<volVectorField> U(movingPhaseModels_[0]*movingPhaseModels_[0].U());

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        U.ref() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].U();
    }

    if (stationaryPhaseModels_.empty())
    {
        return U;
    }
    else
    {
        return U/sumAlphaMoving();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phaseInterfaceKey& key) const
{
    if (interfaceSurfaceTensionModels_.found(key))
    {
        return interfaceSurfaceTensionModels_[key]->sigma();
    }
    else
    {
        return volScalarField::New
        (
            interfaceSurfaceTensionModel::typeName + ":sigma",
            mesh_,
            dimensionedScalar(interfaceSurfaceTensionModel::dimSigma, 0)
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::phaseSystem::sigma(const phaseInterfaceKey& key, const label patchi) const
{
    if (interfaceSurfaceTensionModels_.found(key))
    {
        return interfaceSurfaceTensionModels_[key]->sigma(patchi);
    }
    else
    {
        return tmp<scalarField>
        (
            new scalarField(mesh_.boundary()[patchi].size(), 0)
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        volScalarField::New
        (
            "nearInterface",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    forAll(phases(), phasei)
    {
        tnearInt.ref() = max
        (
            tnearInt(),
            pos0(phases()[phasei] - 0.01)*pos0(0.99 - phases()[phasei])
        );
    }

    return tnearInt;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::dmdtf
(
    const phaseInterfaceKey& key
) const
{
    return volScalarField::New
    (
        IOobject::groupName("dmdtf", phaseInterface(*this, key).name()),
        mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    );
}


Foam::PtrList<Foam::volScalarField> Foam::phaseSystem::dmdts() const
{
    return PtrList<volScalarField>(phaseModels_.size());
}


Foam::PtrList<Foam::volScalarField> Foam::phaseSystem::d2mdtdps() const
{
    return PtrList<volScalarField>(phaseModels_.size());
}


bool Foam::phaseSystem::incompressible() const
{
    forAll(phaseModels_, phasei)
    {
        if (!phaseModels_[phasei].incompressible())
        {
            return false;
        }
    }

    return true;
}


bool Foam::phaseSystem::implicitPhasePressure(const phaseModel& phase) const
{
    return false;
}


bool Foam::phaseSystem::implicitPhasePressure() const
{
    return false;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::surfaceTension
(
    const phaseModel& phase1
) const
{
    tmp<surfaceScalarField> tSurfaceTension
    (
        surfaceScalarField::New
        (
            "surfaceTension",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0)
        )
    );

    forAll(phases(), phasej)
    {
        const phaseModel& phase2 = phases()[phasej];

        if (&phase2 != &phase1)
        {
            const phaseInterface interface(phase1, phase2);

            if (cAlphas_.found(interface))
            {
                tSurfaceTension.ref() +=
                    fvc::interpolate(sigma(interface)*K(phase1, phase2))
                   *(
                        fvc::interpolate(phase2)*fvc::snGrad(phase1)
                      - fvc::interpolate(phase1)*fvc::snGrad(phase2)
                    );
            }
        }
    }

    return tSurfaceTension;
}


void Foam::phaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}


void Foam::phaseSystem::correctContinuityError()
{
    const PtrList<volScalarField> dmdts = this->dmdts();

    forAll(movingPhaseModels_, movingPhasei)
    {
        phaseModel& phase = movingPhaseModels_[movingPhasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.rho();

        volScalarField source
        (
            volScalarField::New
            (
                IOobject::groupName("source", phase.name()),
                mesh_,
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        if (fvModels().addsSupToField(rho.name()))
        {
            source += fvModels().source(alpha, rho)&rho;
        }

        if (dmdts.set(phase.index()))
        {
            source += dmdts[phase.index()];
        }

        phase.correctContinuityError(source);
    }
}


void Foam::phaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctKinematics();

        // Update the pressure time-derivative if required
        if (!updateDpdt && phaseModels_[phasei].thermo().dpdt())
        {
            dpdt_ = fvc::ddt(phaseModels_[phasei].fluidThermo().p());
            updateDpdt = true;
        }
    }
}


void Foam::phaseSystem::correctThermo()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermo();
    }
}


void Foam::phaseSystem::correctReactions()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctReactions();
    }
}


void Foam::phaseSystem::correctSpecies()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctSpecies();
    }
}


void Foam::phaseSystem::predictMomentumTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].predictMomentumTransport();
    }
}


void Foam::phaseSystem::predictThermophysicalTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].predictThermophysicalTransport();
    }
}


void Foam::phaseSystem::correctMomentumTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctMomentumTransport();
    }
}


void Foam::phaseSystem::correctThermophysicalTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermophysicalTransport();
    }
}


void Foam::phaseSystem::meshUpdate()
{
    if (mesh_.changing())
    {
        MRF_.update();

        // forAll(phaseModels_, phasei)
        // {
        //     phaseModels_[phasei].meshUpdate();
        // }
    }
}


void Foam::phaseSystem::correctBoundaryFlux()
{
    forAll(movingPhases(), movingPhasei)
    {
        phaseModel& phase = movingPhases()[movingPhasei];

        tmp<volVectorField> tU(phase.U());
        const volVectorField::Boundary& UBf = tU().boundaryField();

        FieldField<fvsPatchField, scalar> phiRelBf
        (
            MRF_.relative(mesh_.Sf().boundaryField() & UBf)
        );

        surfaceScalarField::Boundary& phiBf = phase.phiRef().boundaryFieldRef();

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                isA<fixedValueFvsPatchScalarField>(phiBf[patchi])
             && !isA<movingWallVelocityFvPatchVectorField>(UBf[patchi])
            )
            {
                phiBf[patchi] == phiRelBf[patchi];
            }
        }
    }
}


void Foam::phaseSystem::correctPhi
(
    const volScalarField& p_rgh,
    const autoPtr<volScalarField>& divU,
    const pressureReference& pressureReference,
    nonOrthogonalSolutionControl& pimple
)
{
    forAll(movingPhases(), movingPhasei)
    {
        phaseModel& phase = movingPhases()[movingPhasei];

        volVectorField::Boundary& Ubf = phase.URef().boundaryFieldRef();
        surfaceVectorField::Boundary& UfBf = phase.UfRef().boundaryFieldRef();

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].initEvaluate();
            }
        }

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].evaluate();
                UfBf[patchi] = Ubf[patchi];
            }
        }
    }

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    correctBoundaryFlux();

    phi_ = Zero;
    PtrList<surfaceScalarField> alphafs(phaseModels_.size());
    forAll(movingPhases(), movingPhasei)
    {
        phaseModel& phase = movingPhases()[movingPhasei];
        const label phasei = phase.index();
        const volScalarField& alpha = phase;

        alphafs.set(phasei, fvc::interpolate(alpha).ptr());

        // Calculate absolute flux
        // from the mapped surface velocity
        phi_ += alphafs[phasei]*(mesh_.Sf() & phase.UfRef());
    }

    if (incompressible())
    {
        fv::correctPhi
        (
            phi_,
            movingPhases()[0].U(),
            p_rgh,
            autoPtr<volScalarField>(),
            divU,
            pressureReference,
            pimple
        );
    }
    else
    {
        volScalarField psi
        (
            volScalarField::New
            (
                "psi",
                mesh_,
                dimensionedScalar(dimless/dimPressure, 0)
            )
        );

        forAll(phases(), phasei)
        {
            phaseModel& phase = phases()[phasei];
            const volScalarField& alpha = phase;

            psi += alpha*phase.fluidThermo().psi()/phase.rho();
        }

        fv::correctPhi
        (
            phi_,
            p_rgh,
            psi,
            autoPtr<volScalarField>(),
            divU(),
            pimple
        );
    }

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi_, movingPhases()[0].U());

    setMixturePhi(alphafs, phi_);
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        forAll(phaseModels_, phasei)
        {
            readOK &= phaseModels_[phasei].read();
        }

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField> Foam::byDt(const volScalarField& vf)
{
    if (fv::localEulerDdt::enabled(vf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaT(vf.mesh())*vf;
    }
    else
    {
        return vf/vf.mesh().time().deltaT();
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::byDt(const surfaceScalarField& sf)
{
    if (fv::localEulerDdt::enabled(sf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaTf(sf.mesh())*sf;
    }
    else
    {
        return sf/sf.mesh().time().deltaT();
    }
}


// ************************************************************************* //
