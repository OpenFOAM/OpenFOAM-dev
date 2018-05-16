/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "populationBalanceModel.H"
#include "coalescenceModel.H"
#include "breakupModel.H"
#include "binaryBreakupModel.H"
#include "driftModel.H"
#include "nucleationModel.H"
#include "phaseSystem.H"
#include "fvmDdt.H"
#include "fvcDdt.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void
Foam::diameterModels::populationBalanceModel::registerVelocityAndSizeGroups()
{
    forAll(fluid_.phases(), phasei)
    {
        if (isA<velocityGroup>(fluid_.phases()[phasei].dPtr()()))
        {
            const velocityGroup& velGroup =
                refCast<const velocityGroup>(fluid_.phases()[phasei].dPtr()());

            if (velGroup.popBalName() == this->name())
            {
                velocityGroups_.append(&const_cast<velocityGroup&>(velGroup));

                forAll(velGroup.sizeGroups(), i)
                {
                    this->add
                    (
                        &const_cast<sizeGroup&>(velGroup.sizeGroups()[i])
                    );
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::add(sizeGroup* group)
{
    if
    (
        sizeGroups_.size() != 0
        &&
        group->x().value() <= sizeGroups_.last()->x().value()
    )
    {
        FatalErrorIn
        (
            "populationBalance::add"
            "(sizeGroup* group)"
        )   << "Size groups must be entered according to their representative"
            << " size"
            << endl
            << exit(FatalError);
    }

    sizeGroups_.append(group);

    // Grid generation over property space
    if (sizeGroups_.size() == 1)
    {
        // Set the first sizeGroup boundary to the representative value of
        // the first sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last()->x()
            )
        );

        // Set the last sizeGroup boundary to the representative size of the
        // last sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last()->x()
            )
        );
    }
    else
    {
        // Overwrite the next-to-last sizeGroup boundary to lie halfway between
        // the last two sizeGroups
        v_.last() =
            0.5
           *(
                sizeGroups_[sizeGroups_.size()-2]->x()
              + sizeGroups_.last()->x()
            );

        // Set the last sizeGroup boundary to the representative size of the
        // last sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last()->x()
            )
        );
    }

    delta_.append(new PtrList<dimensionedScalar>());

    Su_.append
    (
        new volScalarField
        (
            IOobject
            (
                "Su",
                fluid_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "Su",
                inv(dimTime),
                0.0
            )
        )
    );

    SuSp_.append
    (
        new volScalarField
        (
            IOobject
            (
                "SuSp",
                fluid_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "SuSp",
                inv(dimTime),
                0.0
            )
        )
    );
}


void Foam::diameterModels::populationBalanceModel::createPhasePairs()
{
    forAll(velocityGroups_, i)
    {
        const phaseModel& phasei = velocityGroups_[i]->phase();

        forAll(velocityGroups_, j)
        {
            const phaseModel& phasej = velocityGroups_[j]->phase();

            if (&phasei != &phasej)
            {
                const phasePairKey key
                (
                    phasei.name(),
                    phasej.name(),
                    false
                );

                if (!phasePairs_.found(key))
                {
                    phasePairs_.insert
                    (
                        key,
                        autoPtr<phasePair>
                        (
                            new phasePair
                            (
                                phasei,
                                phasej
                            )
                        )
                    );
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::correct()
{
    calcDeltas();

    forAll(velocityGroups_, v)
    {
        velocityGroups_[v]->preSolve();
    }

    forAll(coalescence_, model)
    {
        coalescence_[model].correct();
    }

    forAll(breakup_, model)
    {
        breakup_[model].correct();

        breakup_[model].dsdPtr()().correct();
    }

    forAll(binaryBreakup_, model)
    {
        binaryBreakup_[model].correct();
    }

    forAll(drift_, model)
    {
        drift_[model].correct();
    }

    forAll(nucleation_, model)
    {
        nucleation_[model].correct();
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByCoalescence
(
    const label j,
    const label k
)
{
    const sizeGroup& fj = *sizeGroups_[j];
    const sizeGroup& fk = *sizeGroups_[k];

    dimensionedScalar Gamma;
    dimensionedScalar v = fj.x() + fk.x();

    for (label i = j; i < sizeGroups_.size(); i++)
    {
        // Calculate fraction for intra-interval events
        Gamma = gamma(i, v);

        if (Gamma.value() == 0) continue;

        const sizeGroup& fi = *sizeGroups_[i];

        // Avoid double counting of events
        if (j == k)
        {
            Sui_ =
                0.5*fi.x()*coalescenceRate_()*Gamma
               *fj*fj.phase()/fj.x()
               *fk*fk.phase()/fk.x();
        }
        else
        {
            Sui_ =
                fi.x()*coalescenceRate_()*Gamma
               *fj*fj.phase()/fj.x()
               *fk*fk.phase()/fk.x();
        }

        Su_[i] += Sui_;

        dimensionedScalar ratio = fj.x()/fi.x();

        const phasePairKey pairij
        (
            fi.phase().name(),
            fj.phase().name()
        );

        if (pDmdt_.found(pairij))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
            );

            pDmdt_[pairij]->ref() += dmdtSign*ratio*Sui_*fi.phase().rho();
        }

        const phasePairKey pairik
        (
            fi.phase().name(),
            fk.phase().name()
        );

        if (pDmdt_.found(pairik))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairik).key(), pairik)
            );

            pDmdt_[pairik]->ref() += dmdtSign*(1 - ratio)*Sui_*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::
deathByCoalescence
(
    const label i,
    const label j
)
{
    const sizeGroup& fi = *sizeGroups_[i];
    const sizeGroup& fj = *sizeGroups_[j];

    SuSp_[i] += coalescenceRate_()*fi.phase()*fj*fj.phase()/fj.x();

    if (i != j)
    {
        SuSp_[j] += coalescenceRate_()*fj.phase()*fi*fi.phase()/fi.x();
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByBreakup
(
    const label k,
    const label model
)
{
    const sizeGroup& fk = *sizeGroups_[k];

    for (label i = 0; i <= k; i++)
    {
        const sizeGroup& fi = *sizeGroups_[i];

        Sui_ =
            fi.x()*breakupRate_()*breakup_[model].dsdPtr()().nik(i, k)
           *fk*fk.phase()/fk.x();

        Su_[i] += Sui_;

        const phasePairKey pair
        (
            fi.phase().name(),
            fk.phase().name()
        );

        if (pDmdt_.found(pair))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pair).key(), pair)
            );

            pDmdt_[pair]->ref() += dmdtSign*Sui_*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBreakup(const label i)
{
    const sizeGroup& fi = *sizeGroups_[i];

    SuSp_[i] += breakupRate_()*fi.phase();
}


void Foam::diameterModels::populationBalanceModel::calcDeltas()
{
    forAll(sizeGroups_, i)
    {
        if (delta_[i].empty())
        {
            for (label j = 0; j <= sizeGroups_.size() - 1; j++)
            {
                delta_[i].append
                (
                    new dimensionedScalar
                    (
                        "delta",
                        dimVolume,
                        v_[i+1].value() - v_[i].value()
                    )
                );

                if
                (
                    v_[i].value() < 0.5*sizeGroups_[j]->x().value()
                 &&
                    0.5*sizeGroups_[j]->x().value() < v_[i+1].value()
                )
                {
                    delta_[i][j] =  mag(0.5*sizeGroups_[j]->x() - v_[i]);
                }
                else if (0.5*sizeGroups_[j]->x().value() <= v_[i].value())
                {
                    delta_[i][j] *= 0.0;
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByBinaryBreakup
(
    const label i,
    const label j
)
{
    const sizeGroup& fj = *sizeGroups_[j];
    const sizeGroup& fi = *sizeGroups_[i];

    Sui_ = fi.x()*binaryBreakupRate_()*delta_[i][j]*fj*fj.phase()/fj.x();

    Su_[i] += Sui_;

    const phasePairKey pairij
    (
        fi.phase().name(),
        fj.phase().name()
    );

    if (pDmdt_.found(pairij))
    {
        const scalar dmdtSign
        (
            Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
        );

        pDmdt_[pairij]->ref() += dmdtSign*Sui_*fi.phase().rho();
    }

    dimensionedScalar Gamma;
    dimensionedScalar v = fj.x() - fi.x();

    for (label k = 0; k <= j; k++)
    {
        // Calculate fraction for intra-interval events
        Gamma = gamma(k, v);

        if (Gamma.value() == 0) continue;

        const sizeGroup& fk = *sizeGroups_[k];

        volScalarField& Suk = Sui_;

        Suk =
            sizeGroups_[k]->x()*binaryBreakupRate_()*delta_[i][j]*Gamma
           *fj*fj.phase()/fj.x();

        Su_[k] += Suk;

        const phasePairKey pairkj
        (
            fk.phase().name(),
            fj.phase().name()
        );

        if (pDmdt_.found(pairkj))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare
                (
                    pDmdt_.find(pairkj).key(),
                    pairkj
                )
            );

            pDmdt_[pairkj]->ref() += dmdtSign*Suk*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::
deathByBinaryBreakup
(
    const label j,
    const label i
)
{
    const volScalarField& alphai = sizeGroups_[i]->phase();

    SuSp_[i] += alphai*binaryBreakupRate_()*delta_[j][i];
}


void Foam::diameterModels::populationBalanceModel::drift(const label i)
{
    const sizeGroup& fi = *sizeGroups_[i];

    if (i == 0)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+1]->x()/sizeGroups_[i]->x();
    }
    else if (i == sizeGroups_.size() - 1)
    {
        rx_() = neg(driftRate_())*sizeGroups_[i-1]->x()/sizeGroups_[i]->x();
    }
    else
    {
        rx_() =
            pos(driftRate_())*sizeGroups_[i+1]->x()/sizeGroups_[i]->x()
          + neg(driftRate_())*sizeGroups_[i-1]->x()/sizeGroups_[i]->x();
    }

    SuSp_[i] +=
        (neg(1 - rx_()) + neg(1 - rx_()/(1 - rx_())))*driftRate_()
       *fi.phase()/((rx_() - 1)*sizeGroups_[i]->x());

    rx_() *= 0.0;
    rdx_() *= 0.0;

    if (i < sizeGroups_.size() - 2)
    {
        rx_() += pos(driftRate_())*sizeGroups_[i+2]->x()/sizeGroups_[i+1]->x();

        rdx_() +=
            pos(driftRate_())
           *(sizeGroups_[i+2]->x() - sizeGroups_[i+1]->x())
           /(sizeGroups_[i+1]->x() - sizeGroups_[i]->x());
    }
    else if (i == sizeGroups_.size() - 2)
    {
        rx_() += pos(driftRate_())*sizeGroups_[i+1]->x()/sizeGroups_[i]->x();

        rdx_() +=
            pos(driftRate_())
           *(sizeGroups_[i+1]->x() - sizeGroups_[i]->x())
           /(sizeGroups_[i]->x() - sizeGroups_[i-1]->x());
    }

    if (i == 1)
    {
        rx_() +=
            neg(driftRate_())*sizeGroups_[i-1]->x()
           /sizeGroups_[i]->x();

        rdx_() +=
            neg(driftRate_())
           *(sizeGroups_[i]->x() - sizeGroups_[i-1]->x())
           /(sizeGroups_[i+1]->x() - sizeGroups_[i]->x());
    }
    else if (i > 1)
    {
        rx_() += neg(driftRate_())*sizeGroups_[i-2]->x()/sizeGroups_[i-1]->x();

        rdx_() +=
            neg(driftRate_())
           *(sizeGroups_[i-1]->x() - sizeGroups_[i-2]->x())
           /(sizeGroups_[i]->x() - sizeGroups_[i-1]->x());
    }

    if (i != sizeGroups_.size() - 1)
    {
        const sizeGroup& fj = *sizeGroups_[i+1];
        volScalarField& Suj = Sui_;

        Suj =
            pos(driftRate_())*driftRate_()*rdx_()
           *fi*fi.phase()/fi.x()
           /(rx_() - 1);

        Su_[i+1] += Suj;

        const phasePairKey pairij
        (
            fi.phase().name(),
            fj.phase().name()
        );

        if (pDmdt_.found(pairij))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
            );

            pDmdt_[pairij]->ref() -= dmdtSign*Suj*fi.phase().rho();
        }
    }

    if (i != 0)
    {
        const sizeGroup& fh = *sizeGroups_[i-1];
        volScalarField& Suh = Sui_;

        Suh =
            neg(driftRate_())*driftRate_()*rdx_()
           *fi*fi.phase()/fi.x()
           /(rx_() - 1);

        Su_[i-1] += Suh;

        const phasePairKey pairih
        (
            fi.phase().name(),
            fh.phase().name()
        );

        if (pDmdt_.found(pairih))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairih).key(), pairih)
            );

            pDmdt_[pairih]->ref() -= dmdtSign*Suh*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::nucleation(const label i)
{
    dimensionedScalar volume("volume", dimVolume, 1.0);

    Su_[i] += sizeGroups_[i]->x()*nucleationRate_();
}


void Foam::diameterModels::populationBalanceModel::sources()
{
    forAll(sizeGroups_, i)
    {
        Su_[i] *= 0.0;
        SuSp_[i] *= 0.0;
    }

    forAllConstIter
    (
        phasePairTable,
        phasePairs(),
        phasePairIter
    )
    {
        pDmdt_(phasePairIter())->ref() *= 0.0;
    }

    // Since the calculation of the rates is computationally expensive,
    // they are calculated once for each sizeGroup pair and inserted into source
    // terms as required
    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = *sizeGroups_[i];

        if (coalescence_.size() != 0)
        {
            for (label j = 0; j <= i; j++)
            {
                const sizeGroup& fj = *sizeGroups_[j];

                if (fi.x() + fj.x() > sizeGroups_.last()->x()) break;

                coalescenceRate_() *= 0.0;

                forAll(coalescence_, model)
                {
                    coalescence_[model].addToCoalescenceRate
                    (
                        coalescenceRate_(),
                        i,
                        j
                    );
                }

                birthByCoalescence(i, j);

                deathByCoalescence(i, j);
            }
        }

        if (breakup_.size() != 0)
        {
            forAll(breakup_, model)
            {
                breakup_[model].setBreakupRate(breakupRate_(), i);

                birthByBreakup(i, model);

                deathByBreakup(i);
            }
        }

        if (binaryBreakup_.size() != 0)
        {
            label j = 0;

            while (delta_[j][i].value() != 0.0)
            {
                binaryBreakupRate_() *= 0.0;

                forAll(binaryBreakup_, model)
                {
                    binaryBreakup_[model].addToBinaryBreakupRate
                    (
                        binaryBreakupRate_(),
                        j,
                        i
                    );
                }

                birthByBinaryBreakup(j, i);

                deathByBinaryBreakup(j, i);

                j++;
            }
        }

        if (drift_.size() != 0)
        {
            driftRate_() *= 0.0;

            forAll(drift_, model)
            {
                drift_[model].addToDriftRate(driftRate_(), i);
            }

            drift(i);
        }

        if (nucleation_.size() != 0)
        {
            nucleationRate_() *= 0.0;

            forAll(nucleation_, model)
            {
                nucleation_[model].addToNucleationRate(nucleationRate_(), i);
            }

            nucleation(i);
        }
    }
}


void Foam::diameterModels::populationBalanceModel::dmdt()
{
    forAll(velocityGroups_, v)
    {
        velocityGroup& VelocityGroup = *velocityGroups_[v];

        velocityGroups_[v]->dmdt() *= 0.0;

        forAll(sizeGroups_, i)
        {
            if (&sizeGroups_[i]->phase() == &VelocityGroup.phase())
            {
                sizeGroup& fi = *sizeGroups_[i];

                VelocityGroup.dmdt() += fi.phase().rho()*(Su_[i] - SuSp_[i]*fi);
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::calcAlphas()
{
    alphas_ *= 0.0;

    forAll(velocityGroups_, v)
    {
        alphas_ += velocityGroups_[v]->phase();
    }
}


void Foam::diameterModels::populationBalanceModel::calcVelocity()
{
    U_ *= 0.0;

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v]->phase();

        U_ += phase*phase.U()/max(alphas_, phase.residualAlpha());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceModel::populationBalanceModel
(
    const phaseSystem& fluid,
    const word& name,
    HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>& pDmdt
)
:
    regIOobject
    (
        IOobject
        (
            name,
            fluid.mesh().time().constant(),
            fluid.mesh()
        )
    ),
    fluid_(fluid),
    pDmdt_(pDmdt),
    mesh_(fluid.mesh()),
    name_(name),
    dict_
    (
        fluid.subDict("populationBalanceCoeffs").subDict(name_)
    ),
    pimple_(mesh_.lookupObject<pimpleControl>("solutionControl")),
    continuousPhase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict_.lookup("continuousPhase"))
        )
    ),
    velocityGroups_(),
    sizeGroups_(),
    v_(),
    delta_(),
    Su_(),
    SuSp_(),
    Sui_
    (
        IOobject
        (
            "Sui",
            fluid.time().timeName(),
            fluid.mesh()
        ),
        mesh_,
        dimensionedScalar("Sui", inv(dimTime), Zero)
    ),
    coalescence_
    (
        dict_.lookup("coalescenceModels"),
        coalescenceModel::iNew(*this)
    ),
    coalescenceRate_(),
    breakup_
    (
        dict_.lookup("breakupModels"),
        breakupModel::iNew(*this)
    ),
    breakupRate_(),
    binaryBreakup_
    (
        dict_.lookup("binaryBreakupModels"),
        binaryBreakupModel::iNew(*this)
    ),
    binaryBreakupRate_(),
    drift_
    (
        dict_.lookup("driftModels"),
        driftModel::iNew(*this)
    ),
    driftRate_(),
    rx_(),
    rdx_(),
    nucleation_
    (
        dict_.lookup("nucleationModels"),
        nucleationModel::iNew(*this)
    ),
    nucleationRate_(),
    alphas_
    (
        IOobject
        (
            IOobject::groupName("alpha", this->name()),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            IOobject::groupName("alpha", this->name()),
            dimless,
            Zero
        )
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedVector
        (
            IOobject::groupName("U", this->name()),
            dimVelocity,
            Zero
        )
    )
{
    this->registerVelocityAndSizeGroups();

    this->createPhasePairs();

    if (sizeGroups_.size() < 3)
    {
        FatalErrorInFunction
                    << "The populationBalance " << name_
                    << " requires a minimum number of three sizeGroups to be"
                    << " specified."
                    << exit(FatalError);
    }


    if (coalescence_.size() != 0)
    {
        coalescenceRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                     "coalescenceRate",
                     fluid.time().timeName(),
                     fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("coalescenceRate", dimVolume/dimTime, Zero)
            )
        );
    }

    if (breakup_.size() != 0)
    {
        breakupRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "breakupRate",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("breakupRate", inv(dimTime), Zero)
            )
        );
    }

    if (binaryBreakup_.size() != 0)
    {
        binaryBreakupRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "binaryBreakupRate",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar
                (
                    "binaryBreakupRate",
                    inv(dimVolume*dimTime),
                    Zero
                )
            )
        );
    }

    if (drift_.size() != 0)
    {
        driftRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "driftRate",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("driftRate", dimVolume/dimTime, Zero)
            )
        );

        rx_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "r",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("r", dimless, Zero)
            )
        );

        rdx_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "r",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("r", dimless, Zero)
            )
        );
    }

    if (nucleation_.size() != 0)
    {
        nucleationRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "nucleationRate",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar
                (
                    "nucleationRate",
                    inv(dimTime*dimVolume),
                    Zero
                )
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceModel::~populationBalanceModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::populationBalanceModel>
Foam::diameterModels::populationBalanceModel::clone() const
{
    notImplemented("populationBalance::clone() const");
    return autoPtr<populationBalanceModel>(nullptr);
}


bool Foam::diameterModels::populationBalanceModel::writeData(Ostream& os) const
{
    return os.good();
}


const Foam::dimensionedScalar
Foam::diameterModels::populationBalanceModel::gamma
(
    const label i,
    const dimensionedScalar& v
) const
{
    dimensionedScalar lowerBoundary(v);
    dimensionedScalar upperBoundary(v);
    const dimensionedScalar& xi = sizeGroups_[i]->x();

    if (i == 0)
    {
       lowerBoundary = xi;
    }
    else
    {
       lowerBoundary = sizeGroups_[i-1]->x();
    }

    if (i == sizeGroups_.size() - 1)
    {
        upperBoundary = xi;
    }
    else
    {
        upperBoundary = sizeGroups_[i+1]->x();
    }

    if (v < lowerBoundary || v > upperBoundary)
    {
        return 0.0;
    }
    else if (v.value() <= xi.value())
    {
        return (v - lowerBoundary)/(xi - lowerBoundary);
    }
    else
    {
        return (upperBoundary - v)/(upperBoundary - xi);
    }
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceModel::sigmaWithContinuousPhase
(
    const phaseModel& dispersedPhase
) const
{
    const phasePairKey key
    (
        dispersedPhase.name(),
        continuousPhase_.name()
    );

    return fluid_.sigma(key);
}


const Foam::phaseCompressibleTurbulenceModel&
Foam::diameterModels::populationBalanceModel::continuousTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                continuousPhase_.name()
            )
        );
}


void Foam::diameterModels::populationBalanceModel::solve()
{
    const dictionary& solutionControls = mesh_.solverDict(name_);
    bool solveOnFinalIterOnly
        (
            solutionControls.lookupOrDefault<bool>
            (
                "solveOnFinalIterOnly",
                false
            )
        );

    calcAlphas();
    calcVelocity();

    if (!solveOnFinalIterOnly || pimple_.finalIter())
    {
        label nCorr(readLabel(solutionControls.lookup("nCorr")));
        scalar tolerance
            (
                readScalar(solutionControls.lookup("tolerance"))
            );

        if (nCorr > 0)
        {
            correct();
        }

        int iCorr = 0;
        scalar maxInitialResidual = 1;

        while (++iCorr <= nCorr && maxInitialResidual > tolerance)
        {
            Info<< "populationBalance "
                << this->name()
                << ": Iteration "
                << iCorr
                << endl;

            sources();

            dmdt();

            maxInitialResidual = 0;

            forAll(sizeGroups_, i)
            {
                sizeGroup& fi = *sizeGroups_[i];
                const phaseModel& phase = fi.phase();
                const volScalarField& alpha = phase;
                const dimensionedScalar& residualAlpha = phase.residualAlpha();
                const volScalarField& rho = phase.thermo().rho();

                fvScalarMatrix sizeGroupEqn
                (
                    fvm::ddt(alpha, rho, fi)
                  + fi.VelocityGroup().mvConvection()->fvmDiv
                    (
                        phase.alphaRhoPhi(),
                        fi
                    )
                  - fvm::Sp
                    (
                        fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                      - fi.VelocityGroup().dmdt(),
                        fi
                    )
                  ==
                    Su_[i]*rho
                  - fvm::SuSp(SuSp_[i]*rho, fi)
                  + fvc::ddt(residualAlpha*rho, fi)
                  - fvm::ddt(residualAlpha*rho, fi)
                );

                sizeGroupEqn.relax();

                maxInitialResidual = max
                (
                    sizeGroupEqn.solve().initialResidual(),
                    maxInitialResidual
                );
            }
        }

        if (nCorr > 0)
        {
            forAll(velocityGroups_, i)
            {
                velocityGroups_[i]->postSolve();
            }
        }

        volScalarField fAlpha0
        (
            *sizeGroups_.first()*sizeGroups_.first()->phase()
        );

        volScalarField fAlphaN
        (
            *sizeGroups_.last()*sizeGroups_.last()->phase()
        );

        Info<< this->name() << " sizeGroup phase fraction first, last = "
            << fAlpha0.weightedAverage(this->mesh().V()).value()
            << ' ' << fAlphaN.weightedAverage(this->mesh().V()).value()
            << endl;
    }
}

// ************************************************************************* //
