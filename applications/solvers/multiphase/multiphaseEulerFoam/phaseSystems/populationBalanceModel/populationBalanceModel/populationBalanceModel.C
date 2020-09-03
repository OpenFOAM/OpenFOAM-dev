/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
#include "surfaceTensionModel.H"
#include "fvm.H"
#include "fvcDdt.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "shapeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(populationBalanceModel, 0);
}
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::diameterModels::populationBalanceModel::registerVelocityGroups()
{
    forAll(fluid_.phases(), phasei)
    {
        if (isA<velocityGroup>(fluid_.phases()[phasei].dPtr()()))
        {
            const velocityGroup& velGroup =
                refCast<const velocityGroup>(fluid_.phases()[phasei].dPtr()());

            if (velGroup.popBalName() == this->name())
            {
                velocityGroups_.resize(velocityGroups_.size() + 1);

                velocityGroups_.set
                (
                    velocityGroups_.size() - 1,
                    &const_cast<velocityGroup&>(velGroup)
                );

                forAll(velGroup.sizeGroups(), i)
                {
                    this->registerSizeGroups
                    (
                        const_cast<sizeGroup&>(velGroup.sizeGroups()[i])
                    );
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::registerSizeGroups
(
    sizeGroup& group
)
{
    if
    (
        sizeGroups_.size() != 0
        &&
        group.x().value() <= sizeGroups_.last().x().value()
    )
    {
        FatalErrorInFunction
            << "Size groups must be entered according to their representative"
            << " size"
            << exit(FatalError);
    }

    sizeGroups_.resize(sizeGroups_.size() + 1);
    sizeGroups_.set(sizeGroups_.size() - 1, &group);

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
                sizeGroups_.last().x()
            )
        );

        // Set the last sizeGroup boundary to the representative size of the
        // last sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last().x()
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
                sizeGroups_[sizeGroups_.size()-2].x()
              + sizeGroups_.last().x()
            );

        // Set the last sizeGroup boundary to the representative size of the
        // last sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last().x()
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
            dimensionedScalar(inv(dimTime), 0)
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
            dimensionedScalar(inv(dimTime), 0)
        )
    );
}


void Foam::diameterModels::populationBalanceModel::createPhasePairs()
{
    forAll(velocityGroups_, i)
    {
        const phaseModel& phasei = velocityGroups_[i].phase();

        forAll(velocityGroups_, j)
        {
            const phaseModel& phasej = velocityGroups_[j].phase();

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


void Foam::diameterModels::populationBalanceModel::precompute()
{
    calcDeltas();

    forAll(coalescence_, model)
    {
        coalescence_[model].precompute();
    }

    forAll(breakup_, model)
    {
        breakup_[model].precompute();

        breakup_[model].dsdPtr()->precompute();
    }

    forAll(binaryBreakup_, model)
    {
        binaryBreakup_[model].precompute();
    }

    forAll(drift_, model)
    {
        drift_[model].precompute();
    }

    forAll(nucleation_, model)
    {
        nucleation_[model].precompute();
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByCoalescence
(
    const label j,
    const label k
)
{
    const sizeGroup& fj = sizeGroups_[j];
    const sizeGroup& fk = sizeGroups_[k];

    dimensionedScalar Eta;
    dimensionedScalar v = fj.x() + fk.x();

    for (label i = j; i < sizeGroups_.size(); i++)
    {
        // Calculate fraction for intra-interval events
        Eta = eta(i, v);

        if (Eta.value() == 0) continue;

        const sizeGroup& fi = sizeGroups_[i];

        // Avoid double counting of events
        if (j == k)
        {
            Sui_ =
                0.5*fi.x()*coalescenceRate_()*Eta
               *fj*fj.phase()/fj.x()
               *fk*fk.phase()/fk.x();
        }
        else
        {
            Sui_ =
                fi.x()*coalescenceRate_()*Eta
               *fj*fj.phase()/fj.x()
               *fk*fk.phase()/fk.x();
        }

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

            *pDmdt_[pairij] += dmdtSign*fj.x()/v*Sui_*fi.phase().rho();
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

            *pDmdt_[pairik] += dmdtSign*fk.x()/v*Sui_*fi.phase().rho();
        }

        sizeGroups_[i].shapeModelPtr()->addCoalescence(Sui_, fj, fk);
    }
}


void Foam::diameterModels::populationBalanceModel::
deathByCoalescence
(
    const label i,
    const label j
)
{
    const sizeGroup& fi = sizeGroups_[i];
    const sizeGroup& fj = sizeGroups_[j];

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
    const sizeGroup& fk = sizeGroups_[k];

    for (label i = 0; i <= k; i++)
    {
        const sizeGroup& fi = sizeGroups_[i];

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

            *pDmdt_[pair] += dmdtSign*Sui_*fi.phase().rho();
        }

        sizeGroups_[i].shapeModelPtr()->addBreakup(Sui_, fk);
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBreakup(const label i)
{
    const sizeGroup& fi = sizeGroups_[i];

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
                    v_[i].value() < 0.5*sizeGroups_[j].x().value()
                 &&
                    0.5*sizeGroups_[j].x().value() < v_[i+1].value()
                )
                {
                    delta_[i][j] =  mag(0.5*sizeGroups_[j].x() - v_[i]);
                }
                else if (0.5*sizeGroups_[j].x().value() <= v_[i].value())
                {
                    delta_[i][j].value() = 0;
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
    const sizeGroup& fj = sizeGroups_[j];
    const sizeGroup& fi = sizeGroups_[i];

    Sui_ = fi.x()*binaryBreakupRate_()*delta_[i][j]*fj*fj.phase()/fj.x();

    Su_[i] += Sui_;

    sizeGroups_[i].shapeModelPtr()->addBreakup(Sui_, fj);

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

        *pDmdt_[pairij] += dmdtSign*Sui_*fi.phase().rho();
    }

    dimensionedScalar Eta;
    dimensionedScalar v = fj.x() - fi.x();

    for (label k = 0; k <= j; k++)
    {
        // Calculate fraction for intra-interval events
        Eta = eta(k, v);

        if (Eta.value() == 0) continue;

        const sizeGroup& fk = sizeGroups_[k];

        volScalarField& Suk = Sui_;

        Suk =
            sizeGroups_[k].x()*binaryBreakupRate_()*delta_[i][j]*Eta
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

            *pDmdt_[pairkj] += dmdtSign*Suk*fi.phase().rho();
        }

        sizeGroups_[i].shapeModelPtr()->addBreakup(Suk, fj);
    }
}


void Foam::diameterModels::populationBalanceModel::
deathByBinaryBreakup
(
    const label j,
    const label i
)
{
    const volScalarField& alphai = sizeGroups_[i].phase();

    SuSp_[i] += alphai*binaryBreakupRate_()*delta_[j][i];
}


void Foam::diameterModels::populationBalanceModel::drift
(
    const label i,
    driftModel& model
)
{
    model.addToDriftRate(driftRate_(), i);

    const sizeGroup& fp = sizeGroups_[i];

    if (i == 0)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+1].x()/sizeGroups_[i].x();
    }
    else if (i == sizeGroups_.size() - 1)
    {
        rx_() = neg(driftRate_())*sizeGroups_[i-1].x()/sizeGroups_[i].x();
    }
    else
    {
        rx_() =
            pos(driftRate_())*sizeGroups_[i+1].x()/sizeGroups_[i].x()
          + neg(driftRate_())*sizeGroups_[i-1].x()/sizeGroups_[i].x();
    }

    SuSp_[i] +=
        (neg(1 - rx_()) + neg(rx_() - rx_()/(1 - rx_())))*driftRate_()
       *fp.phase()/((rx_() - 1)*fp.x());

    rx_() = Zero;
    rdx_() = Zero;

    if (i == sizeGroups_.size() - 2)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+1].x()/sizeGroups_[i].x();

        rdx_() =
            pos(driftRate_())
           *(sizeGroups_[i+1].x() - sizeGroups_[i].x())
           /(sizeGroups_[i].x() - sizeGroups_[i-1].x());
    }
    else if (i < sizeGroups_.size() - 2)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+2].x()/sizeGroups_[i+1].x();

        rdx_() =
            pos(driftRate_())
           *(sizeGroups_[i+2].x() - sizeGroups_[i+1].x())
           /(sizeGroups_[i+1].x() - sizeGroups_[i].x());
    }

    if (i == 1)
    {
        rx_() += neg(driftRate_())*sizeGroups_[i-1].x()/sizeGroups_[i].x();

        rdx_() +=
            neg(driftRate_())
           *(sizeGroups_[i].x() - sizeGroups_[i-1].x())
           /(sizeGroups_[i+1].x() - sizeGroups_[i].x());
    }
    else if (i > 1)
    {
        rx_() += neg(driftRate_())*sizeGroups_[i-2].x()/sizeGroups_[i-1].x();

        rdx_() +=
            neg(driftRate_())
           *(sizeGroups_[i-1].x() - sizeGroups_[i-2].x())
           /(sizeGroups_[i].x() - sizeGroups_[i-1].x());
    }

    if (i != sizeGroups_.size() - 1)
    {
        const sizeGroup& fe = sizeGroups_[i+1];
        volScalarField& Sue = Sui_;

        Sue =
            pos(driftRate_())*driftRate_()*rdx_()
           *fp*fp.phase()/fp.x()
           /(rx_() - 1);

        Su_[i+1] += Sue;

        const phasePairKey pairij
        (
            fp.phase().name(),
            fe.phase().name()
        );

        if (pDmdt_.found(pairij))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
            );

            *pDmdt_[pairij] -= dmdtSign*Sue*fp.phase().rho();
        }

        sizeGroups_[i+1].shapeModelPtr()->addDrift(Sue, fp, model);
    }

    if (i != 0)
    {
        const sizeGroup& fw = sizeGroups_[i-1];
        volScalarField& Suw = Sui_;

        Suw =
            neg(driftRate_())*driftRate_()*rdx_()
           *fp*fp.phase()/fp.x()
           /(rx_() - 1);

        Su_[i-1] += Suw;

        const phasePairKey pairih
        (
            fp.phase().name(),
            fw.phase().name()
        );

        if (pDmdt_.found(pairih))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairih).key(), pairih)
            );

            *pDmdt_[pairih] -= dmdtSign*Suw*fp.phase().rho();
        }

        sizeGroups_[i-1].shapeModelPtr()->addDrift(Suw, fp, model);
    }
}


void Foam::diameterModels::populationBalanceModel::
nucleation
(
    const label i,
    nucleationModel& model
)
{
    const sizeGroup& fi = sizeGroups_[i];

    model.addToNucleationRate(nucleationRate_(), i);

    Sui_ = sizeGroups_[i].x()*nucleationRate_();

    Su_[i] += Sui_;

    sizeGroups_[i].shapeModelPtr()->addNucleation(Sui_, fi, model);
}


void Foam::diameterModels::populationBalanceModel::sources()
{
    forAll(sizeGroups_, i)
    {
        sizeGroups_[i].shapeModelPtr()->reset();
        Su_[i] = Zero;
        SuSp_[i] = Zero;
    }

    forAllConstIter
    (
        phasePairTable,
        phasePairs(),
        phasePairIter
    )
    {
        *pDmdt_(phasePairIter()) = Zero;
    }

    forAll(sizeGroups_, i)
    {
        if (coalescence_.size() != 0)
        {
            for (label j = 0; j <= i; j++)
            {
                coalescenceRate_() = Zero;

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

            while (delta_[j][i].value() != 0)
            {
                binaryBreakupRate_() = Zero;

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
            forAll(drift_, model)
            {
                driftRate_() = Zero;
                drift(i, drift_[model]);
            }
        }

        if (nucleation_.size() != 0)
        {
            forAll(nucleation_, model)
            {
                nucleationRate_() = Zero;
                nucleation(i, nucleation_[model]);
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::calcAlphas()
{
    alphas_() = Zero;

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v].phase();

        alphas_() += max(phase, phase.residualAlpha());
    }
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceModel::calcDsm()
{
    tmp<volScalarField> tInvDsm
    (
        volScalarField::New
        (
            "invDsm",
            mesh_,
            dimensionedScalar(inv(dimLength), Zero)
        )
    );

    volScalarField& invDsm = tInvDsm.ref();

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v].phase();

        invDsm += max(phase, phase.residualAlpha())/(phase.d()*alphas_());
    }

    return 1.0/tInvDsm;
}


void Foam::diameterModels::populationBalanceModel::calcVelocity()
{
    U_() = Zero;

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v].phase();

        U_() += max(phase, phase.residualAlpha())*phase.U()/alphas_();
    }
}

bool Foam::diameterModels::populationBalanceModel::updateSources()
{
    const bool result = sourceUpdateCounter_ % sourceUpdateInterval() == 0;

    ++ sourceUpdateCounter_;

    return result;
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
            fluid.time().constant(),
            fluid.mesh()
        )
    ),
    fluid_(fluid),
    pDmdt_(pDmdt),
    mesh_(fluid_.mesh()),
    name_(name),
    dict_
    (
        fluid_.subDict("populationBalanceCoeffs").subDict(name_)
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
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), Zero)
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
    alphas_(),
    dsm_(),
    U_(),
    sourceUpdateCounter_(0)
{
    this->registerVelocityGroups();

    this->createPhasePairs();

    if (sizeGroups_.size() < 3)
    {
        FatalErrorInFunction
            << "The populationBalance " << name_
            << " requires a minimum number of three sizeGroups to be specified."
            << exit(FatalError);
    }


    if (coalescence_.size() != 0)
    {
        coalescenceRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                     "coalescenceRate",
                     mesh_.time().timeName(),
                     mesh_
                ),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
            )
        );
    }

    if (breakup_.size() != 0)
    {
        breakupRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "breakupRate",
                    fluid_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(inv(dimTime), Zero)
            )
        );
    }

    if (binaryBreakup_.size() != 0)
    {
        binaryBreakupRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "binaryBreakupRate",
                    fluid_.time().timeName(),
                    mesh_
                ),
                mesh_,
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
        driftRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "driftRate",
                    fluid_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
            )
        );

        rx_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "r",
                    fluid_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );

        rdx_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "r",
                    fluid_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }

    if (nucleation_.size() != 0)
    {
        nucleationRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "nucleationRate",
                    fluid_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar
                (
                    "nucleationRate",
                    inv(dimTime*dimVolume),
                    Zero
                )
            )
        );
    }

    if (velocityGroups_.size() > 1)
    {
        alphas_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", this->name()),
                    fluid_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );

        dsm_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("d", this->name()),
                    fluid_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimLength, Zero)
            )
        );

        U_.set
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("U", this->name()),
                    fluid_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimVelocity, Zero)
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
Foam::diameterModels::populationBalanceModel::eta
(
    const label i,
    const dimensionedScalar& v
) const
{
    const dimensionedScalar& x0 = sizeGroups_[0].x();
    const dimensionedScalar& xi = sizeGroups_[i].x();
    const dimensionedScalar& xm = sizeGroups_.last().x();
    dimensionedScalar lowerBoundary(x0);
    dimensionedScalar upperBoundary(xm);

    if (i != 0) lowerBoundary = sizeGroups_[i-1].x();

    if (i != sizeGroups_.size() - 1) upperBoundary = sizeGroups_[i+1].x();

    if ((i == 0 && v < x0) || (i == sizeGroups_.size() - 1 && v > xm))
    {
        return v/xi;
    }
    else if (v < lowerBoundary || v > upperBoundary)
    {
        return 0;
    }
    else if (v.value() == xi.value())
    {
        return 1.0;
    }
    else if (v > xi)
    {
        return (upperBoundary - v)/(upperBoundary - xi);
    }
    else
    {
        return (v - lowerBoundary)/(xi - lowerBoundary);
    }
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceModel::sigmaWithContinuousPhase
(
    const phaseModel& dispersedPhase
) const
{
    return
        fluid_.lookupSubModel<surfaceTensionModel>
        (
            phasePair(dispersedPhase, continuousPhase_)
        ).sigma();
}


const Foam::phaseCompressibleMomentumTransportModel&
Foam::diameterModels::populationBalanceModel::continuousTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressibleMomentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                continuousPhase_.name()
            )
        );
}


void Foam::diameterModels::populationBalanceModel::solve()
{
    const dictionary& solutionControls = mesh_.solverDict(name_);
    bool solveOnFinalIterOnly =
        solutionControls.lookupOrDefault<bool>("solveOnFinalIterOnly", false);

    if (!solveOnFinalIterOnly || pimple_.finalPimpleIter())
    {
        const label nCorr = this->nCorr();
        const scalar tolerance =
            solutionControls.lookup<scalar>("tolerance");

        if (nCorr > 0)
        {
            precompute();
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

            if (updateSources())
            {
                sources();
            }

            maxInitialResidual = 0;

            forAll(sizeGroups_, i)
            {
                sizeGroup& fi = sizeGroups_[i];
                const phaseModel& phase = fi.phase();
                const volScalarField& alpha = phase;
                const dimensionedScalar& residualAlpha = phase.residualAlpha();
                const volScalarField& rho = phase.thermo().rho();

                fvScalarMatrix sizeGroupEqn
                (
                    fvm::ddt(alpha, fi)
                  + fvm::div(phase.alphaPhi(), fi)
                ==
                    Su_[i]
                  - fvm::SuSp(SuSp_[i], fi)
                  + fluid_.fvOptions()(alpha, rho, fi)/rho
                  + fvc::ddt(residualAlpha, fi)
                  - fvm::ddt(residualAlpha, fi)
                );

                sizeGroupEqn.relax();
                fluid_.fvOptions().constrain(sizeGroupEqn);

                maxInitialResidual = max
                (
                    sizeGroupEqn.solve().initialResidual(),
                    maxInitialResidual
                );

                fluid_.fvOptions().correct(fi);
            }
        }

        volScalarField fAlpha0
        (
            sizeGroups_.first()*sizeGroups_.first().phase()
        );

        volScalarField fAlphaN
        (
            sizeGroups_.last()*sizeGroups_.last().phase()
        );

        Info<< this->name() << " sizeGroup phase fraction first, last = "
            << fAlpha0.weightedAverage(this->mesh().V()).value()
            << ' ' << fAlphaN.weightedAverage(this->mesh().V()).value()
            << endl;
    }
}


void Foam::diameterModels::populationBalanceModel::correct()
{
    if (velocityGroups_.size() > 1)
    {
        calcAlphas();
        dsm_() = calcDsm();
        calcVelocity();
    }
}


// ************************************************************************* //
