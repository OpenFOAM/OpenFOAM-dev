/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "ThermalPhaseChangePhaseSystem.H"
#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fvcVolumeIntegrate.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
addMomentumTransfer(phaseSystem::momentumTransferTable& eqns) const
{
    // Source term due to mass trasfer
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const volVectorField& U1(pair.phase1().U());
        const volVectorField& U2(pair.phase2().U());

        const volScalarField dmdt(this->iDmdt(pair) + this->wDmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[pair.phase1().name()] += dmdt21*U2 - fvm::Sp(dmdt21, U1);
        *eqns[pair.phase2().name()] -= dmdt12*U1 - fvm::Sp(dmdt12, U2);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
ThermalPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    volatile_(this->template lookupOrDefault<word>("volatile", "none")),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    ),
    massTransfer_(this->lookup("massTransfer"))
{

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        // Initially assume no mass transfer
        iDmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

        // Initially assume no mass transfer
        wDmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("wDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

        // Initially assume no mass transfer
        wMDotL_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("wMDotL", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
~ThermalPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::iDmdt
(
    const phasePairKey& key
) const
{
    const scalar dmdtSign(Pair<word>::compare(iDmdt_.find(key).key(), key));
    return dmdtSign**iDmdt_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::iDmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> tiDmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("iDmdt", phase.name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        if (pair.contains(phase))
        {
            tiDmdt.ref() += this->iDmdt
            (
                phasePairKey(phase.name(), pair.otherPhase(phase).name())
            );
        }
    }

    return tiDmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::wDmdt
(
    const phasePairKey& key
) const
{
    const scalar dmdtSign(Pair<word>::compare(wDmdt_.find(key).key(), key));

    return dmdtSign**wDmdt_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::wDmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> twDmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("wDmdt", phase.name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        if (pair.contains(phase))
        {
            twDmdt.ref() += this->wDmdt
            (
                phasePairKey(phase.name(), pair.otherPhase(phase).name())
            );
        }
    }

    return twDmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> tDmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Dmdt", phase.name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        if (pair.contains(phase))
        {
            tDmdt.ref() += this->dmdt
            (
                phasePairKey(phase.name(), pair.otherPhase(phase).name())
            );
        }
    }

    return tDmdt;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    addMomentumTransfer(eqnsPtr());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::momentumTransferf() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransferf());

    addMomentumTransfer(eqnsPtr());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    // Source term due to mass transfer
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        if
        (
            this->heatTransferModels_.found(phasePairIter.key())
        )
        {
            const phasePair& pair(phasePairIter());

            if (pair.ordered())
            {
                continue;
            }

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

            const volScalarField& he1(phase1.thermo().he());
            const volScalarField& he2(phase2.thermo().he());

            const volScalarField& K1(phase1.K());
            const volScalarField& K2(phase2.K());

            const volScalarField dmdt(this->dmdt(pair));
            const volScalarField dmdt21(posPart(dmdt));
            const volScalarField dmdt12(negPart(dmdt));
            const volScalarField& Tf(*this->Tf_[pair]);

            *eqns[phase1.name()] +=
                dmdt21*(phase1.thermo().he(phase1.thermo().p(), Tf))
              - fvm::Sp(dmdt21, he1)
              + dmdt21*(K2 - K1);

            *eqns[phase2.name()] -=
                dmdt12*(phase2.thermo().he(phase2.thermo().p(), Tf))
              - fvm::Sp(dmdt12, he2)
              + dmdt12*(K1 - K2);

            *eqns[phase1.name()] += negPart(*this->wMDotL_[pair]);
            *eqns[phase2.name()] -= posPart(*this->wMDotL_[pair]);
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr =
        BasePhaseSystem::massTransfer();

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {

            if (Yi[i].member() != volatile_)
            {
                continue;
            }

            const word name
            (
                IOobject::groupName(volatile_, phase.name())
            );

            const word otherName
            (
                IOobject::groupName(volatile_, otherPhase.name())
            );

            const volScalarField dmdt(this->iDmdt(pair) + this->wDmdt(pair));
            const volScalarField dmdt12(posPart(dmdt));
            const volScalarField dmdt21(negPart(dmdt));

            *eqns[name] += fvm::Sp(dmdt21, eqns[name]->psi()) - dmdt21;
            *eqns[otherName] +=
                dmdt12 - fvm::Sp(dmdt12, eqns[otherName]->psi());
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::correctThermo()
{
    typedef compressible::alphatPhaseChangeWallFunctionFvPatchScalarField
        alphatPhaseChangeWallFunction;

    phaseSystem::correctThermo();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        if
        (
            this->heatTransferModels_.found(phasePairIter.key())
        )
        {
            const phasePair& pair(phasePairIter());

            if (pair.ordered())
            {
                continue;
            }

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

            const volScalarField& T1(phase1.thermo().T());
            const volScalarField& T2(phase2.thermo().T());

            const volScalarField& he1(phase1.thermo().he());
            const volScalarField& he2(phase2.thermo().he());

            volScalarField& iDmdt(*this->iDmdt_[pair]);
            volScalarField& Tf(*this->Tf_[pair]);

            volScalarField hef1(phase1.thermo().he(phase1.thermo().p(), Tf));
            volScalarField hef2(phase2.thermo().he(phase2.thermo().p(), Tf));

            volScalarField L
            (
                (neg0(iDmdt)*hef2 + pos(iDmdt)*he2)
              - (pos0(iDmdt)*hef1 + neg(iDmdt)*he1)
            );

            volScalarField iDmdtNew(iDmdt);

            if (massTransfer_ )
            {
                volScalarField H1
                (
                    this->heatTransferModels_[pair][pair.first()]->K(0)
                );

                volScalarField H2
                (
                    this->heatTransferModels_[pair][pair.second()]->K(0)
                );

                Tf = saturationModel_->Tsat(phase1.thermo().p());

                iDmdtNew =
                    (H1*(Tf - T1) + H2*(Tf - T2))/L;
            }
            else
            {
                iDmdtNew == dimensionedScalar("0", iDmdt.dimensions(), 0);
            }

            volScalarField H1
            (
                this->heatTransferModels_[pair][pair.first()]->K()
            );

            volScalarField H2
            (
                this->heatTransferModels_[pair][pair.second()]->K()
            );

            // Limit the H[12] to avoid /0
            H1.max(small);
            H2.max(small);

            Tf = (H1*T1 + H2*T2 + iDmdtNew*L)/(H1 + H2);

            Info<< "Tf." << pair.name()
                << ": min = " << min(Tf.primitiveField())
                << ", mean = " << average(Tf.primitiveField())
                << ", max = " << max(Tf.primitiveField())
                << endl;

            scalar iDmdtRelax(this->mesh().fieldRelaxationFactor("iDmdt"));
            iDmdt = (1 - iDmdtRelax)*iDmdt + iDmdtRelax*iDmdtNew;

            if (massTransfer_ )
            {
                Info<< "iDmdt." << pair.name()
                    << ": min = " << min(iDmdt.primitiveField())
                    << ", mean = " << average(iDmdt.primitiveField())
                    << ", max = " << max(iDmdt.primitiveField())
                    << ", integral = " << fvc::domainIntegrate(iDmdt).value()
                    << endl;
            }

            volScalarField& wDmdt(*this->wDmdt_[pair]);
            volScalarField& wMDotL(*this->wMDotL_[pair]);
            wDmdt *= 0;
            wMDotL *= 0;

            bool wallBoilingActive = false;

            forAllConstIter(phasePair, pair, iter)
            {
                const phaseModel& phase = iter();
                const phaseModel& otherPhase = iter.otherPhase();

                if
                (
                    phase.mesh().foundObject<volScalarField>
                    (
                        "alphat." +  phase.name()
                    )
                )
                {
                    const volScalarField& alphat =
                        phase.mesh().lookupObject<volScalarField>
                        (
                            "alphat." +  phase.name()
                        );

                    const fvPatchList& patches = this->mesh().boundary();
                    forAll(patches, patchi)
                    {
                        const fvPatch& currPatch = patches[patchi];

                        if
                        (
                            isA<alphatPhaseChangeWallFunction>
                            (
                                alphat.boundaryField()[patchi]
                            )
                        )
                        {
                            const alphatPhaseChangeWallFunction& PCpatch =
                                refCast<const alphatPhaseChangeWallFunction>
                                (
                                    alphat.boundaryField()[patchi]
                                );

                            phasePairKey key(phase.name(), otherPhase.name());

                            if (PCpatch.activePhasePair(key))
                            {
                                wallBoilingActive = true;

                                const scalarField& patchDmdt =
                                    PCpatch.dmdt(key);
                                const scalarField& patchMDotL =
                                    PCpatch.mDotL(key);

                                const scalar sign
                                (
                                    Pair<word>::compare(pair, key)
                                );

                                forAll(patchDmdt, facei)
                                {
                                    const label faceCelli =
                                        currPatch.faceCells()[facei];
                                    wDmdt[faceCelli] -= sign*patchDmdt[facei];
                                    wMDotL[faceCelli] -= sign*patchMDotL[facei];
                                }
                            }
                        }
                    }
                }
            }

            if (wallBoilingActive)
            {
                Info<< "wDmdt." << pair.name()
                    << ": min = " << min(wDmdt.primitiveField())
                    << ", mean = " << average(wDmdt.primitiveField())
                    << ", max = " << max(wDmdt.primitiveField())
                    << ", integral = " << fvc::domainIntegrate(wDmdt).value()
                    << endl;
            }
        }
    }
}


template<class BasePhaseSystem>
bool Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
