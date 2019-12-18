/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::totalDmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tTotalDmdtf = phaseSystem::dmdtf(key);
    volScalarField& totalDmdtf = tTotalDmdtf.ref();

    const scalar tTotalDmdtfSign =
        Pair<word>::compare(this->phasePairs_[key], key);

    if (dmdtfs_.found(key))
    {
        totalDmdtf += *dmdtfs_[key];
    }

    if (nDmdtfs_.found(key))
    {
        totalDmdtf += *nDmdtfs_[key];
    }

    return tTotalDmdtfSign*tTotalDmdtf;
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
    phaseChange_(this->lookup("phaseChange"))
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
        dmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:dmdtf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        // Initially assume no mass transfer
        nDmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:nucleation:dmdtf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        // Initially assume no mass transfer
        nDmdtLfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:nucleation:dmdtLf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
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
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    return BasePhaseSystem::dmdtf(key) + this->totalDmdtf(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[dmdtfIter.key()];
        const volScalarField& dmdtf = *dmdtfIter();

        addField(pair.phase1(), "dmdt", dmdtf, dmdts);
        addField(pair.phase2(), "dmdt", - dmdtf, dmdts);
    }

    forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[nDmdtfIter.key()];
        const volScalarField& nDmdtf = *nDmdtfIter();

        addField(pair.phase1(), "dmdt", nDmdtf, dmdts);
        addField(pair.phase2(), "dmdt", - nDmdtf, dmdts);
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);
    this->addDmdtUfs(nDmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);
    this->addDmdtUfs(nDmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        typename BasePhaseSystem::heatTransferModelTable,
        this->heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& Tf(*this->Tf_[pair]);

        const volScalarField H1(heatTransferModelIter().first()->K());
        const volScalarField H2(heatTransferModelIter().second()->K());
        const volScalarField HEff(H1*H2/(H1 + H2));

        *eqns[phase1.name()] +=
            H1*(Tf - phase1.thermo().T())
          - HEff*(phase2.thermo().T() - phase1.thermo().T());

        *eqns[phase2.name()] +=
            H2*(Tf - phase2.thermo().T())
          - HEff*(phase1.thermo().T() - phase2.thermo().T());
    }

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

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        const volScalarField dmdtf(this->totalDmdtf(pair));
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        *eqns[phase1.name()] += - fvm::Sp(dmdtf21, he1) + dmdtf21*(K2 - K1);

        *eqns[phase2.name()] -= - fvm::Sp(dmdtf12, he2) + dmdtf12*(K1 - K2);

        if (this->heatTransferModels_.found(phasePairIter.key()))
        {
            const volScalarField& Tf(*this->Tf_[pair]);

            *eqns[phase1.name()] +=
                dmdtf21*phase1.thermo().he(phase1.thermo().p(), Tf);

            *eqns[phase2.name()] -=
                dmdtf12*phase2.thermo().he(phase2.thermo().p(), Tf);
        }
        else
        {
            *eqns[phase1.name()] += dmdtf21*he2;

            *eqns[phase2.name()] -= dmdtf12*he1;
        }

        if (this->nDmdtLfs_.found(phasePairIter.key()))
        {
            *eqns[phase1.name()] += negPart(*this->nDmdtLfs_[pair]);
            *eqns[phase2.name()] -= posPart(*this->nDmdtLfs_[pair]);

            if
            (
                phase1.thermo().he().member() == "e"
             || phase2.thermo().he().member() == "e"
            )
            {
                if (phase1.thermo().he().member() == "e")
                {
                    *eqns[phase1.name()] +=
                        phase1.thermo().p()*dmdtf/phase1.thermo().rho();
                }

                if (phase2.thermo().he().member() == "e")
                {
                    *eqns[phase2.name()] -=
                        phase2.thermo().p()*dmdtf/phase2.thermo().rho();
                }
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

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

        if (!phase.pure() || !otherPhase.pure())
        {
            FatalErrorInFunction
                << "ThermalPhaseChangePhaseSystem does not currently support "
                << "multiComponent phase models."
                << exit(FatalError);
        }

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

            // Note that the phase YiEqn does not contain a continuity error
            // term, so these additions represent the entire mass transfer

            const volScalarField dmdtf(this->totalDmdtf(pair));

            *eqns[name] += dmdtf;
            *eqns[otherName] -= dmdtf;
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::correctInterfaceThermo()
{
    typedef compressible::alphatPhaseChangeWallFunctionFvPatchScalarField
        alphatPhaseChangeWallFunction;

    forAllConstIter
    (
        typename BasePhaseSystem::heatTransferModelTable,
        this->heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& T1(phase1.thermo().T());
        const volScalarField& T2(phase2.thermo().T());

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField& p(phase1.thermo().p());

        volScalarField& dmdtf(*this->dmdtfs_[pair]);
        volScalarField& Tf(*this->Tf_[pair]);

        const volScalarField Tsat(saturationModel_->Tsat(phase1.thermo().p()));

        volScalarField hf1
        (
            he1.member() == "e"
          ? phase1.thermo().he(p, Tsat) + p/phase1.rho()
          : phase1.thermo().he(p, Tsat)
        );
        volScalarField hf2
        (
            he2.member() == "e"
          ? phase2.thermo().he(p, Tsat) + p/phase2.rho()
          : phase2.thermo().he(p, Tsat)
        );

        volScalarField h1
        (
            he1.member() == "e"
          ? he1 + p/phase1.rho()
          : tmp<volScalarField>(he1)
        );

        volScalarField h2
        (
            he2.member() == "e"
          ? he2 + p/phase2.rho()
          : tmp<volScalarField>(he2)
        );

        volScalarField L
        (
            (neg0(dmdtf)*hf2 + pos(dmdtf)*h2)
          - (pos0(dmdtf)*hf1 + neg(dmdtf)*h1)
        );

        volScalarField dmdtfNew(dmdtf);

        if (phaseChange_)
        {
            volScalarField H1(heatTransferModelIter().first()->K(0));
            volScalarField H2(heatTransferModelIter().second()->K(0));

            dmdtfNew = (H1*(Tsat - T1) + H2*(Tsat - T2))/L;
        }
        else
        {
            dmdtfNew == dimensionedScalar(dmdtf.dimensions(), 0);
        }

        volScalarField H1(heatTransferModelIter().first()->K());
        volScalarField H2(heatTransferModelIter().second()->K());

        // Limit the H[12] to avoid /0
        H1.max(small);
        H2.max(small);

        Tf = (H1*T1 + H2*T2 + dmdtfNew*L)/(H1 + H2);

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.primitiveField())
            << ", mean = " << average(Tf.primitiveField())
            << ", max = " << max(Tf.primitiveField())
            << endl;

        scalar dmdtfRelax =
            this->mesh().fieldRelaxationFactor(dmdtf.member());
        dmdtf = (1 - dmdtfRelax)*dmdtf + dmdtfRelax*dmdtfNew;

        if (phaseChange_)
        {
            Info<< dmdtf.name()
                << ": min = " << min(dmdtf.primitiveField())
                << ", mean = " << average(dmdtf.primitiveField())
                << ", max = " << max(dmdtf.primitiveField())
                << ", integral = " << fvc::domainIntegrate(dmdtf).value()
                << endl;
        }

        volScalarField& nDmdtf(*this->nDmdtfs_[pair]);
        volScalarField& nDmdtLf(*this->nDmdtLfs_[pair]);
        nDmdtf = Zero;
        nDmdtLf = Zero;

        bool wallBoilingActive = false;

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            if
            (
                phase.mesh().foundObject<volScalarField>
                (
                    IOobject::groupName("alphat", phase.name())
                )
            )
            {
                const volScalarField& alphat =
                    phase.mesh().lookupObject<volScalarField>
                    (
                        IOobject::groupName("alphat", phase.name())
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

                            const scalarField& patchDmdtf = PCpatch.dmdtf(key);
                            const scalarField& patchDmdtLf =
                                PCpatch.dmdtLf(key);

                            const scalar sign
                            (
                                Pair<word>::compare(pair, key)
                            );

                            forAll(patchDmdtf, facei)
                            {
                                const label faceCelli =
                                    currPatch.faceCells()[facei];
                                nDmdtf[faceCelli] -= sign*patchDmdtf[facei];
                                nDmdtLf[faceCelli] -= sign*patchDmdtLf[facei];
                            }
                        }
                    }
                }
            }
        }

        if (wallBoilingActive)
        {
            Info<< nDmdtf.name()
                << ": min = " << min(nDmdtf.primitiveField())
                << ", mean = " << average(nDmdtf.primitiveField())
                << ", max = " << max(nDmdtf.primitiveField())
                << ", integral = " << fvc::domainIntegrate(nDmdtf).value()
                << endl;
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
