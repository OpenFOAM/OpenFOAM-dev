/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
#include "rhoReactionThermo.H"

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
    volatile_(this->template lookupOrDefault<word>("volatile", "none"))
{
    this->generatePairsAndSubModels
    (
        "saturation",
        saturationModels_
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
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::saturation
(
    const phasePairKey& key
) const
{
    return saturationModels_[key];
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

        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();

        const volScalarField& he1(thermo1.he());
        const volScalarField& he2(thermo2.he());

        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        const volScalarField dmdtf(this->totalDmdtf(pair));
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        *eqns[phase1.name()] += - fvm::Sp(dmdtf21, he1) + dmdtf21*(K2 - K1);

        *eqns[phase2.name()] -= - fvm::Sp(dmdtf12, he2) + dmdtf12*(K1 - K2);

        if (this->saturationModels_.found(phasePairIter.key()))
        {
            const volScalarField& Tf(*this->Tf_[pair]);

            if (volatile_ != "none" && isA<rhoReactionThermo>(thermo1))
            {
                const basicSpecieMixture& composition1 =
                    refCast<const rhoReactionThermo>(thermo1).composition();

                *eqns[phase1.name()] +=
                    dmdtf21
                   *composition1.HE
                    (
                        composition1.species()[volatile_],
                        thermo1.p(),
                        Tf
                    );
            }
            else
            {
                *eqns[phase1.name()] +=
                    dmdtf21*thermo1.he(thermo1.p(), Tf);
            }

            if (volatile_ != "none" && isA<rhoReactionThermo>(thermo2))
            {
                const basicSpecieMixture& composition2 =
                    refCast<const rhoReactionThermo>(thermo2).composition();

                *eqns[phase2.name()] -=
                    dmdtf12
                   *composition2.HE
                    (
                        composition2.species()[volatile_],
                        thermo2.p(),
                        Tf
                    );
            }
            else
            {
                *eqns[phase2.name()] -=
                    dmdtf12*thermo2.he(thermo2.p(), Tf);
            }
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
        }

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

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    if (volatile_ == "none")
    {
        return eqnsPtr;
    }

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

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        if (phase1.pure() || phase2.pure())
        {
            FatalErrorInFunction
                << "Volatile specie was given, but at least one phase in pair "
                << pair << " is pure."
                << exit(FatalError);
        }

        const volScalarField& Y1 = phase1.Y(volatile_);
        const volScalarField& Y2 = phase2.Y(volatile_);

        // Note that the phase YiEqn does not contain a continuity error
        // term, so these additions represent the entire mass transfer

        const volScalarField dmdtf(this->totalDmdtf(pair));

        *eqns[Y1.name()] += dmdtf;
        *eqns[Y2.name()] -= dmdtf;
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

        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();

        const volScalarField& T1(thermo1.T());
        const volScalarField& T2(thermo2.T());

        const volScalarField& p(thermo1.p());

        volScalarField& dmdtf(*this->dmdtfs_[pair]);
        volScalarField& Tf(*this->Tf_[pair]);

        volScalarField dmdtfNew(dmdtf);

        if (saturationModels_.found(heatTransferModelIter.key()))
        {
            const phasePairKey& key = heatTransferModelIter.key();
            const volScalarField Tsat(saturation(key).Tsat(thermo1.p()));

            volScalarField ha1(thermo1.ha());
            volScalarField haf1(thermo1.ha(p, Tsat));

            if (volatile_ != "none" && isA<rhoReactionThermo>(thermo1))
            {
                 const basicSpecieMixture& composition1 =
                     refCast<const rhoReactionThermo>(thermo1).composition();

                 ha1 =
                     composition1.Ha
                     (
                         composition1.species()[volatile_],
                         p,
                         T1
                     );

                 haf1 =
                     composition1.Ha
                     (
                         composition1.species()[volatile_],
                         p,
                         Tsat
                     );
             }

             volScalarField ha2(thermo2.ha());
             volScalarField haf2(thermo2.ha(p, Tsat));

             if (volatile_ != "none" && isA<rhoReactionThermo>(thermo2))
             {
                 const basicSpecieMixture& composition2 =
                     refCast<const rhoReactionThermo>(thermo2).composition();

                 ha2 =
                     composition2.Ha
                     (
                         composition2.species()[volatile_],
                         p,
                         T2
                     );

                 haf2 =
                     composition2.Ha
                     (
                         composition2.species()[volatile_],
                         p,
                         Tsat
                     );
             }

            volScalarField L
            (
                (neg0(dmdtf)*haf2 + pos(dmdtf)*ha2)
              - (pos0(dmdtf)*haf1 + neg(dmdtf)*ha1)
            );

            volScalarField H1(heatTransferModelIter().first()->K(0));
            volScalarField H2(heatTransferModelIter().second()->K(0));

            dmdtfNew = (H1*(Tsat - T1) + H2*(Tsat - T2))/L;

            if (volatile_ != "none")
            {
                dmdtfNew *=
                    neg0(dmdtfNew)*phase1.Y(volatile_)
                  + pos(dmdtfNew)*phase2.Y(volatile_);
            }

            H1 = heatTransferModelIter().first()->K();
            H2 = heatTransferModelIter().second()->K();

            // Limit the H[12] to avoid /0
            H1.max(small);
            H2.max(small);

            Tf = (H1*T1 + H2*T2 + dmdtfNew*L)/(H1 + H2);
        }
        else
        {
            dmdtfNew = Zero;
            volScalarField H1(heatTransferModelIter().first()->K());
            volScalarField H2(heatTransferModelIter().second()->K());

            Tf = (H1*T1 + H2*T2)/(H1 + H2);
        }

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.primitiveField())
            << ", mean = " << average(Tf.primitiveField())
            << ", max = " << max(Tf.primitiveField())
            << endl;

        scalar dmdtfRelax =
            this->mesh().fieldRelaxationFactor(dmdtf.member());
        dmdtf = (1 - dmdtfRelax)*dmdtf + dmdtfRelax*dmdtfNew;

        if (saturationModels_.found(heatTransferModelIter.key()))
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
                ) &&
                saturationModels_.found(heatTransferModelIter.key())
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
