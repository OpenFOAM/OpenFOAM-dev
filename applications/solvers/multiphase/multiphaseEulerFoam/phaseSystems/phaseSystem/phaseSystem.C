/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
#include "surfaceTensionModel.H"
#include "aspectRatioModel.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"
#include "localEulerDdtScheme.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "CorrectPhi.H"
#include "fvcMeshPhi.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "unitConversion.H"
#include "dragModel.H"
#include "BlendedInterfacialModel.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "pimpleControl.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
    defineRunTimeSelectionTable(phaseSystem, dictionary);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::calcPhi
(
    const phaseModelList& phaseModels
) const
{
    tmp<surfaceScalarField> tmpPhi
    (
        surfaceScalarField::New
        (
            "phi",
            fvc::interpolate(phaseModels[0])*phaseModels[0].phi()
        )
    );

    for (label phasei=1; phasei<phaseModels.size(); phasei++)
    {
        tmpPhi.ref() +=
            fvc::interpolate(phaseModels[phasei])*phaseModels[phasei].phi();
    }

    return tmpPhi;
}


void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {}

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


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


void Foam::phaseSystem::correctContactAngle
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = phase1.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            alphaContactAngleFvPatchScalarField::thetaPropsTable::
                const_iterator tp =
                    acap.thetaProps()
                   .find(phasePairKey(phase1.name(), phase2.name()));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface "
                    << phasePairKey(phase1.name(), phase2.name())
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == phase1.name());

            scalar theta0 = degToRad(tp().theta0(matched));
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > small)
            {
                scalar thetaA = degToRad(tp().thetaA(matched));
                scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    phase1.U()().boundaryField()[patchi].patchInternalField()
                  - phase1.U()().boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + small);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatPatch & AfHatPatch);

            scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::K
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(phase1, phase2);

    correctContactAngle(phase1, phase2, tnHatfv.ref().boundaryFieldRef());

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
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    referencePhaseName_(lookupOrDefault("referencePhase", word::null)),

    phaseModels_
    (
        lookup("phases"),
        phaseModel::iNew(*this, referencePhaseName_)
    ),

    phi_(calcPhi(phaseModels_)),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure/dimTime, 0)
    ),

    MRF_(mesh_),

    cAlphas_(lookupOrDefault("interfaceCompression", cAlphaTable())),

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
    label multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        movingPhasei += !phase.stationary();
        stationaryPhasei += phase.stationary();
        anisothermalPhasei += !phase.isothermal();
        multiComponentPhasei += !phase.pure();
    }
    movingPhaseModels_.resize(movingPhasei);
    stationaryPhaseModels_.resize(stationaryPhasei);
    anisothermalPhaseModels_.resize(anisothermalPhasei);
    multiComponentPhaseModels_.resize(multiComponentPhasei);

    movingPhasei = 0;
    stationaryPhasei = 0;
    anisothermalPhasei = 0;
    multiComponentPhasei = 0;
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
            multiComponentPhaseModels_.set(multiComponentPhasei++, &phase);
        }
    }

    // Write phi
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    // Blending methods
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().keyword(),
            blendingMethod::New
            (
                iter().keyword(),
                iter().dict(),
                phaseModels_.toc()
            )
        );
    }

    // Sub-models
    generatePairsAndSubModels("surfaceTension", surfaceTensionModels_);
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);

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

    forAll(phases(), phasei)
    {
        const volScalarField& alphai = phases()[phasei];
        mesh_.setFluxRequired(alphai.name());
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
Foam::phaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }
    else
    {
        return volScalarField::New
        (
            aspectRatioModel::typeName + ":E",
            mesh_,
            dimensionedScalar(dimless, 1)
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phasePairKey& key) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma();
    }
    else
    {
        return volScalarField::New
        (
            surfaceTensionModel::typeName + ":sigma",
            mesh_,
            dimensionedScalar(surfaceTensionModel::dimSigma, 0)
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::phaseSystem::sigma(const phasePairKey& key, label patchi) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma(patchi);
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
    const phasePairKey& key
) const
{
    const phasePair pair
    (
        phaseModels_[key.first()],
        phaseModels_[key.second()]
    );

    return volScalarField::New
    (
        IOobject::groupName("dmdtf", pair.name()),
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
            phasePairKey key12(phase1.name(), phase2.name());

            cAlphaTable::const_iterator cAlpha(cAlphas_.find(key12));

            if (cAlpha != cAlphas_.end())
            {
                tSurfaceTension.ref() +=
                    fvc::interpolate(sigma(key12)*K(phase1, phase2))
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
        volScalarField& rho = phase.thermoRef().rho();

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

        updateDpdt = updateDpdt || phaseModels_[phasei].thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.begin()().thermo().p());
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


void Foam::phaseSystem::correctTurbulence()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctTurbulence();
    }
}


void Foam::phaseSystem::correctEnergyTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctEnergyTransport();
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

        const volVectorField::Boundary& UBf = phase.U()().boundaryField();

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
    const tmp<volScalarField>& divU,
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

    {
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
            phi_ += alphafs[phasei]*(mesh_.Sf() & phase.Uf());
        }

        CorrectPhi
        (
            phi_,
            movingPhases()[0].U(),
            p_rgh,
            // surfaceScalarField("rAUf", fvc::interpolate(rAU())),
            dimensionedScalar(dimTime/dimDensity, 1),
            divU(),
            pressureReference,
            pimple
        );

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi_, movingPhases()[0].U());

        setMixturePhi(alphafs, phi_);
    }
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
