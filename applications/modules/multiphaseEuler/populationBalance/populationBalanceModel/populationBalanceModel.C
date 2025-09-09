/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "coalescenceModel.H"
#include "breakupModel.H"
#include "binaryBreakupModel.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "shapeModel.H"
#include "distribution.H"
#include "fvSpecificSource.H"
#include "growthFvScalarFieldSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(populationBalanceModel, 0);
}
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

const Foam::dictionary&
Foam::diameterModels::populationBalanceModel::coeffDict() const
{
    return fluid_.subDict("populationBalanceCoeffs").subDict(name_);
}


void
Foam::diameterModels::populationBalanceModel::precomputeCoalescenceAndBreakup()
{
    forAll(coalescenceModels_, model)
    {
        coalescenceModels_[model].precompute();
    }

    forAll(breakupModels_, model)
    {
        breakupModels_[model].precompute();

        breakupModels_[model].dsdPtr()->precompute();
    }

    forAll(binaryBreakupModels_, model)
    {
        binaryBreakupModels_[model].precompute();
    }
}


void Foam::diameterModels::populationBalanceModel::birthByCoalescence
(
    const label j,
    const label k
)
{
    const sizeGroup& fj = sizeGroups()[j];
    const sizeGroup& fk = sizeGroups()[k];

    const dimensionedScalar v = fj.x() + fk.x();

    for (label i = j; i < sizeGroups().size(); i++)
    {
        const dimensionedScalar Eta = eta(i, v);

        if (Eta.value() == 0) continue;

        const sizeGroup& fi = sizeGroups()[i];

        tmp<volScalarField::Internal> tSui;
        if (j == k)
        {
            tSui =
                0.5*fi.x()/(fj.x()*fk.x())*Eta
               *coalescenceRate_()*fj*fj.phase()*fk*fk.phase();
        }
        else
        {
            tSui =
                fi.x()/(fj.x()*fk.x())*Eta
               *coalescenceRate_()*fj*fj.phase()*fk*fk.phase();
        }
        const volScalarField::Internal& Sui = tSui();

        Su_[i] += Sui;

        const phaseInterface interfaceij(fi.phase(), fj.phase());

        if (dmdtfs_.found(interfaceij))
        {
            const scalar dmdtSign =
                interfaceij.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfaceij] += dmdtSign*fj.x()/v*Sui*fj.phase().rho();
        }

        const phaseInterface interfaceik(fi.phase(), fk.phase());

        if (dmdtfs_.found(interfaceik))
        {
            const scalar dmdtSign =
                interfaceik.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfaceik] += dmdtSign*fk.x()/v*Sui*fk.phase().rho();
        }

        sizeGroups_[i].shape().addCoalescence(Sui, fj, fk);
    }
}


void Foam::diameterModels::populationBalanceModel::deathByCoalescence
(
    const label i,
    const label j
)
{
    const sizeGroup& fi = sizeGroups()[i];
    const sizeGroup& fj = sizeGroups()[j];

    Sp_[i] -= coalescenceRate_()*fi.phase()*fj*fj.phase()/fj.x();

    if (i != j)
    {
        Sp_[j] -= coalescenceRate_()*fj.phase()*fi*fi.phase()/fi.x();
    }
}


void Foam::diameterModels::populationBalanceModel::birthByBreakup
(
    const label k,
    const label model
)
{
    const sizeGroup& fk = sizeGroups()[k];

    for (label i = 0; i <= k; i++)
    {
        const sizeGroup& fi = sizeGroups()[i];

        tmp<volScalarField::Internal> tSui =
            fi.x()*breakupModels_[model].dsdPtr()().nik(i, k)/fk.x()
           *breakupRate_()*fk*fk.phase();
        const volScalarField::Internal& Sui = tSui();

        Su_[i] += Sui;

        const phaseInterface interface(fi.phase(), fk.phase());

        if (dmdtfs_.found(interface))
        {
            const scalar dmdtSign =
                interface.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interface] += dmdtSign*Sui*fk.phase().rho();
        }

        sizeGroups_[i].shape().addBreakup(Sui, fk);
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBreakup(const label i)
{
    Sp_[i] -= breakupRate_()*sizeGroups()[i].phase();
}


void Foam::diameterModels::populationBalanceModel::birthByBinaryBreakup
(
    const label i,
    const label j
)
{
    const sizeGroup& fi = sizeGroups()[i];
    const sizeGroup& fj = sizeGroups()[j];

    const volScalarField::Internal Su(binaryBreakupRate_()*fj*fj.phase());

    {
        tmp<volScalarField::Internal> tSui = fi.x()*delta_[i][j]/fj.x()*Su;
        const volScalarField::Internal& Sui = tSui();

        Su_[i] += Sui;

        sizeGroups_[i].shape().addBreakup(Sui, fj);

        const phaseInterface interfaceij(fi.phase(), fj.phase());

        if (dmdtfs_.found(interfaceij))
        {
            const scalar dmdtSign =
                interfaceij.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfaceij] += dmdtSign*Sui*fj.phase().rho();
        }
    }

    const dimensionedScalar v = fj.x() - fi.x();

    for (label k = 0; k <= j; k ++)
    {
        const dimensionedScalar Eta = eta(k, v);

        if (Eta.value() == 0) continue;

        const sizeGroup& fk = sizeGroups()[k];

        tmp<volScalarField::Internal> tSuk = fk.x()*delta_[i][j]*Eta/fj.x()*Su;
        const volScalarField::Internal& Suk = tSuk();

        Su_[k] += Suk;

        const phaseInterface interfacekj(fk.phase(), fj.phase());

        if (dmdtfs_.found(interfacekj))
        {
            const scalar dmdtSign =
                interfacekj.index(fk.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfacekj] += dmdtSign*Suk*fj.phase().rho();
        }

        sizeGroups_[k].shape().addBreakup(Suk, fj);
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBinaryBreakup
(
    const label j,
    const label i
)
{
    Sp_[i] -= sizeGroups()[i].phase()*binaryBreakupRate_()*delta_[j][i];
}


void
Foam::diameterModels::populationBalanceModel::computeCoalescenceAndBreakup()
{
    forAll(sizeGroups(), i)
    {
        sizeGroups_[i].shape().reset();
    }

    forAll(sizeGroups(), i)
    {
        Su_[i] = Zero;
        Sp_[i] = Zero;
    }

    forAllIter(dmdtfTable, dmdtfs_, dmdtfIter)
    {
        *dmdtfIter() = Zero;
    }

    forAll(coalescencePairs_, coalescencePairi)
    {
        label i = coalescencePairs_[coalescencePairi].first();
        label j = coalescencePairs_[coalescencePairi].second();

        coalescenceRate_() = Zero;

        forAll(coalescenceModels_, model)
        {
            coalescenceModels_[model].addToCoalescenceRate
            (
                coalescenceRate_(),
                i,
                j
            );
        }

        birthByCoalescence(i, j);

        deathByCoalescence(i, j);
    }

    forAll(sizeGroups(), i)
    {
        forAll(breakupModels_, model)
        {
            breakupModels_[model].setBreakupRate(breakupRate_(), i);

            birthByBreakup(i, model);

            deathByBreakup(i);
        }
    }

    forAll(binaryBreakupPairs_, binaryBreakupPairi)
    {
        label i = binaryBreakupPairs_[binaryBreakupPairi].first();
        label j = binaryBreakupPairs_[binaryBreakupPairi].second();

        binaryBreakupRate_() = Zero;

        forAll(binaryBreakupModels_, model)
        {
            binaryBreakupModels_[model].addToBinaryBreakupRate
            (
                binaryBreakupRate_(),
                j,
                i
            );
        }

        birthByBinaryBreakup(j, i);

        deathByBinaryBreakup(j, i);
    }
}


void Foam::diameterModels::populationBalanceModel::precomputeExpansion()
{
    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const velocityGroup& velGrp = *iter();
        const phaseModel& phase = velGrp.phase();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        expansionRates_.set
        (
            phase.index(),
          - (fvc::ddt(alpha, rho)()() + fvc::div(phase.alphaRhoPhi())()())/rho()
          + fvc::ddt(alpha)()() + fvc::div(phase.alphaPhi())()()
        );
    }
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::diameterModels::populationBalanceModel::expansionSus
(
    const label i,
    const UPtrList<const volScalarField>& flds
) const
{
    const sizeGroup& fi = sizeGroups()[i];

    auto fiFld = [&](const label deltai)
    {
        return
            flds.empty()
          ? tmp<volScalarField::Internal>(sizeGroups()[i + deltai])
          : sizeGroups()[i + deltai]()*flds[i + deltai]();
    };

    Pair<tmp<volScalarField::Internal>> tSus;

    if (i != 0)
    {
        const sizeGroup& fiMinus1 = sizeGroups()[i - 1];
        const phaseModel& phaseMinus1 = fiMinus1.phase();

        tSus.first() =
            posPart(expansionRates_[phaseMinus1.index()])
           *fi.x()/(fi.x() - fiMinus1.x())
           *fiFld(-1);
    }

    if (i != sizeGroups().size() - 1)
    {
        const sizeGroup& fiPlus1 = sizeGroups()[i + 1];
        const phaseModel& phasePlus1 = fiPlus1.phase();

        tSus.second() =
          - negPart(expansionRates_[phasePlus1.index()])
           *fi.x()/(fiPlus1.x() - fi.x())
           *fiFld(+1);
    }

    return tSus;
}


void Foam::diameterModels::populationBalanceModel::computeExpansion()
{
    forAllIter(dmdtfTable, expansionDmdtfs_, dmdtfIter)
    {
        *dmdtfIter() = Zero;
    }

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const sizeGroup& fi0 = iter()->sizeGroups().last();

        if (fi0.i() == sizeGroups().size() - 1) continue;

        const sizeGroup& fi1 = sizeGroups()[fi0.i() + 1];

        Pair<tmp<volScalarField::Internal>> tSus0 = expansionSus(fi0.i());
        Pair<tmp<volScalarField::Internal>> tSus1 = expansionSus(fi1.i());

        const phaseInterface interface01(fi0.phase(), fi1.phase());
        const scalar sign = interface01.index(fi0.phase()) == 0 ? -1 : +1;

        *expansionDmdtfs_[interface01] +=
            sign
           *(
              - tSus0.second()*fi0.phase().rho()
              + tSus1.first()*fi1.phase().rho()
            );
    }
}


void Foam::diameterModels::populationBalanceModel::precomputeModelSources()
{
    // nothing to do
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::diameterModels::populationBalanceModel::modelSourceRhoSus
(
    const label i,
    const UPtrList<const volScalarField>& flds
) const
{
    const sizeGroup& fi = sizeGroups()[i];
    const velocityGroup& velGrp = fi.group();

    const volScalarField& fldi = flds.empty() ? fi : flds[fi.i()];

    Pair<tmp<volScalarField::Internal>> tRhoSus;

    forAll(tRhoSus, Sui)
    {
        const sizeGroup& fiPopBalEnd =
            Sui == 0 ? sizeGroups().first() : sizeGroups().last();

        if (fi.i() == fiPopBalEnd.i()) continue;

        const sizeGroup& fiVelGrpEnd =
            Sui == 0 ? velGrp.sizeGroups().first() : velGrp.sizeGroups().last();

        if (fi.i() != fiVelGrpEnd.i()) continue;

        const sizeGroup& fiOther = sizeGroups()[i + (Sui == 0 ? -1 : +1)];

        forAll(fluid_.fvModels(), modeli)
        {
            if (!isA<fvSpecificSource>(fluid_.fvModels()[modeli])) continue;

            const fvSpecificSource& source =
                refCast<const fvSpecificSource>(fluid_.fvModels()[modeli]);

            if
            (
                source.addsSupToField
                (
                    fiOther.phase().volScalarField::name()
                )
             && isA<growthFvScalarFieldSource>
                (
                    fldi.sources()[source.name()]
                )
            )
            {
                const growthFvScalarFieldSource& growthSource =
                    refCast<const growthFvScalarFieldSource>
                    (
                        fldi.sources()[source.name()]
                    );

                const volScalarField::Internal S(source.S(fiOther.name()));

                Pair<tmp<volScalarField::Internal>> sourceCoeffs =
                    growthSource.sourceCoeffs(source);

                tRhoSus[Sui] =
                    Sui == 0
                  ? posPart(S)*sourceCoeffs.first()
                  : negPart(S)*sourceCoeffs.second();
            }
        }
    }

    return tRhoSus;
}


void Foam::diameterModels::populationBalanceModel::computeModelSources()
{
    forAllIter(dmdtfTable, modelSourceDmdtfs_, dmdtfIter)
    {
        *dmdtfIter() = Zero;
    }

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const sizeGroup& fi0 = iter()->sizeGroups().last();

        if (fi0.i() == sizeGroups().size() - 1) continue;

        const sizeGroup& fi1 = sizeGroups()[fi0.i() + 1];

        Pair<tmp<volScalarField::Internal>> tRhoSus0 =
            modelSourceRhoSus(fi0.i());
        Pair<tmp<volScalarField::Internal>> tRhoSus1 =
            modelSourceRhoSus(fi1.i());

        const phaseInterface interface01(fi0.phase(), fi1.phase());
        const scalar sign = interface01.index(fi0.phase()) == 0 ? -1 : +1;

        if (tRhoSus0.second().valid())
        {
            *modelSourceDmdtfs_[interface01] -= sign*tRhoSus0.second();
        }

        if (tRhoSus1.first().valid())
        {
            *modelSourceDmdtfs_[interface01] -= sign*tRhoSus1.first();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::computeDilatationErrors()
{
    PtrList<volScalarField::Internal> modelSourceDmdts(fluid_.phases().size());
    forAllConstIter(dmdtfTable, modelSourceDmdtfs_, dmdtfIter)
    {
        const phaseInterface interface(fluid_, dmdtfIter.key());

        addField(interface.phase1(), "dmdt", *dmdtfIter(), modelSourceDmdts);
        addField(interface.phase2(), "dmdt", - *dmdtfIter(), modelSourceDmdts);
    }

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const velocityGroup& velGrp = *iter();
        const phaseModel& phase = velGrp.phase();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        dilatationErrors_.set
        (
            phase.index(),
            fvc::ddt(alpha)()() + fvc::div(phase.alphaPhi())()()
          - (fluid_.fvModels().source(alpha, rho) & rho)()()/rho()
        );

        forAll(velGrp.sizeGroups(), i)
        {
            const sizeGroup& fi = velGrp.sizeGroups()[i];

            dilatationErrors_[phase.index()] -=
                Su_[fi.i()] + expansionSu(fi.i())
              + (Sp_[fi.i()] + expansionSp(fi.i()))*fi;
        }

        if (modelSourceDmdts.set(phase.index()))
        {
            dilatationErrors_[phase.index()] -=
                modelSourceDmdts[phase.index()]/rho;
        }
    }
}


bool Foam::diameterModels::populationBalanceModel::updateSources()
{
    const bool result = sourceUpdateCounter_ % sourceUpdateInterval() == 0;

    ++ sourceUpdateCounter_;

    return result;
}


Foam::Pair<Foam::dimensionedScalar>
Foam::diameterModels::populationBalanceModel::etaCoeffs0(const label i) const
{
    static const dimensionedScalar z(dimless, scalar(0));

    const dimensionedScalar& xi = sizeGroups()[i].x();

    if (i == 0) return Pair<dimensionedScalar>(z, 1/xi);

    const dimensionedScalar& x0 = sizeGroups()[i - 1].x();

    return Pair<dimensionedScalar>(- x0/(xi - x0), 1/(xi - x0));
}


Foam::Pair<Foam::dimensionedScalar>
Foam::diameterModels::populationBalanceModel::etaCoeffs1(const label i) const
{
    static const dimensionedScalar z(dimless, scalar(0));

    const label n = sizeGroups().size();

    const dimensionedScalar& xi = sizeGroups()[i].x();

    if (i == n - 1) return Pair<dimensionedScalar>(z, 1/xi);

    const dimensionedScalar& x1 = sizeGroups()[i + 1].x();

    return Pair<dimensionedScalar>(x1/(x1 - xi), - 1/(x1 - xi));
}


Foam::Pair<Foam::dimensionedScalar>
Foam::diameterModels::populationBalanceModel::etaVCoeffs0(const label i) const
{
    static const dimensionedScalar o(scalar(1)), zV(dimVolume, scalar(0));

    if (i == 0) return Pair<dimensionedScalar>(o, zV);

    const dimensionedScalar& x0 = 1/sizeGroups()[i - 1].x();
    const dimensionedScalar& xi = 1/sizeGroups()[i].x();

    return Pair<dimensionedScalar>(- x0/(xi - x0), 1/(xi - x0));
}


Foam::Pair<Foam::dimensionedScalar>
Foam::diameterModels::populationBalanceModel::etaVCoeffs1(const label i) const
{
    static const dimensionedScalar o(scalar(1)), zV(dimVolume, scalar(0));

    const label n = sizeGroups().size();

    if (i == n - 1) return Pair<dimensionedScalar>(o, zV);

    const dimensionedScalar& xi = 1/sizeGroups()[i].x();
    const dimensionedScalar& x1 = 1/sizeGroups()[i + 1].x();

    return Pair<dimensionedScalar>(x1/(x1 - xi), - 1/(x1 - xi));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceModel::populationBalanceModel
(
    const phaseSystem& fluid,
    const word& name
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
    mesh_(fluid_.mesh()),
    name_(name),
    continuousPhase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", coeffDict().lookup("continuousPhase"))
        )
    ),
    sizeGroups_(),
    v_(),
    delta_(),
    Su_(),
    Sp_(),
    dmdtfs_(),
    expansionDmdtfs_(),
    modelSourceDmdtfs_(),
    expansionRates_(fluid_.phases().size()),
    dilatationErrors_(fluid_.phases().size()),
    coalescenceModels_
    (
        coeffDict().lookup("coalescenceModels"),
        coalescenceModel::iNew(*this)
    ),
    coalescenceRate_(),
    coalescencePairs_(),
    breakupModels_
    (
        coeffDict().lookup("breakupModels"),
        breakupModel::iNew(*this)
    ),
    breakupRate_(),
    binaryBreakupModels_
    (
        coeffDict().lookup("binaryBreakupModels"),
        binaryBreakupModel::iNew(*this)
    ),
    binaryBreakupRate_(),
    binaryBreakupPairs_(),
    alphas_(),
    dsm_(),
    U_(),
    sourceUpdateCounter_(0)
{
    groups::retrieve(*this, velocityGroupPtrs_, sizeGroups_);

    if (sizeGroups().size() < 3)
    {
        FatalErrorInFunction
            << "The populationBalance " << name_
            << " requires a minimum number of three sizeGroups to be specified."
            << exit(FatalError);
    }

    // Create size-group boundaries
    v_.setSize(sizeGroups().size() + 1);
    v_.set(0, new dimensionedScalar("v", sizeGroups()[0].x()));
    for (label i = 1; i < sizeGroups().size(); ++ i)
    {
        v_.set
        (
            i,
            new dimensionedScalar
            (
                "v",
                (sizeGroups()[i-1].x() + sizeGroups()[i].x())/2
            )
        );
    }
    v_.set(v_.size() - 1, new dimensionedScalar("v", sizeGroups().last().x()));

    // Create section widths if needed
    if (binaryBreakupModels_.size() != 0)
    {
        delta_.setSize(sizeGroups().size());

        forAll(sizeGroups(), i)
        {
            delta_.set(i, new PtrList<dimensionedScalar>());
        }

        forAll(sizeGroups(), i)
        {
            if (!delta_[i].empty()) continue;

            for (label j = 0; j <= sizeGroups().size() - 1; j++)
            {
                const sizeGroup& fj = sizeGroups()[j];

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
                    v_[i].value() < 0.5*fj.x().value()
                 && 0.5*fj.x().value() < v_[i+1].value()
                )
                {
                    delta_[i][j] =  mag(0.5*fj.x() - v_[i]);
                }
                else if (0.5*fj.x().value() <= v_[i].value())
                {
                    delta_[i][j].value() = 0;
                }
            }
        }
    }

    // Create size-group source terms
    Su_.setSize(sizeGroups().size());
    Sp_.setSize(sizeGroups().size());
    forAll(sizeGroups(), i)
    {
        Su_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::groupName("Su" + Foam::name(i), this->name()),
                    fluid_.time().name(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(inv(dimTime), 0)
            )
        );

        Sp_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::groupName("Sp" + Foam::name(i), this->name()),
                    fluid_.time().name(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(inv(dimTime), 0)
            )
        );
    }

    // Create interfacial mass transfer rates
    forAllConstIter
    (
        HashTable<const diameterModels::velocityGroup*>,
        velocityGroupPtrs_,
        iter1
    )
    {
        const diameterModels::velocityGroup& velGrp1 = *iter1();

        forAllConstIter
        (
            HashTable<const diameterModels::velocityGroup*>,
            velocityGroupPtrs_,
            iter2
        )
        {
            const diameterModels::velocityGroup& velGrp2 = *iter2();

            const phaseInterface interface(velGrp1.phase(), velGrp2.phase());

            if (&velGrp1 != &velGrp2 && !dmdtfs_.found(interface))
            {
                dmdtfs_.insert
                (
                    interface,
                    new volScalarField::Internal
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                typedName("dmdtf"),
                                interface.name()
                            ),
                            mesh().time().name(),
                            mesh()
                        ),
                        mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );

                expansionDmdtfs_.insert
                (
                    interface,
                    new volScalarField::Internal
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                typedName("expansionDmdtf"),
                                interface.name()
                            ),
                            mesh().time().name(),
                            mesh()
                        ),
                        mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );

                modelSourceDmdtfs_.insert
                (
                    interface,
                    new volScalarField::Internal
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                typedName("modelSourceDmdtf"),
                                interface.name()
                            ),
                            mesh().time().name(),
                            mesh()
                        ),
                        mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );
            }
        }
    }

    if (coalescenceModels_.size() != 0)
    {
        coalescenceRate_.set
        (
            new volScalarField::Internal
            (
                IOobject
                (
                     IOobject::groupName("coalescenceRate", this->name()),
                     mesh_.time().name(),
                     mesh_
                ),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
            )
        );

        forAll(sizeGroups(), i)
        {
            for (label j = 0; j <= i; j++)
            {
                coalescencePairs_.append(labelPair(i, j));
            }
        }
    }

    if (breakupModels_.size() != 0)
    {
        breakupRate_.set
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::groupName("breakupRate", this->name()),
                    fluid_.time().name(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(inv(dimTime), Zero)
            )
        );
    }

    if (binaryBreakupModels_.size() != 0)
    {
        binaryBreakupRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("binaryBreakupRate", this->name()),
                    fluid_.time().name(),
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

        forAll(sizeGroups(), i)
        {
            label j = 0;

            while (delta_[j][i].value() != 0)
            {
                binaryBreakupPairs_.append(labelPair(i, j));
                j++;
            }
        }
    }

    correct();
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


Foam::dimensionedScalar Foam::diameterModels::populationBalanceModel::eta
(
    const label i,
    const dimensionedScalar& v
) const
{
    const Pair<dimensionedScalar> coeffs0 = etaCoeffs0(i);
    const Pair<dimensionedScalar> coeffs1 = etaCoeffs1(i);

    return
        max
        (
            min
            (
                coeffs0.first() + coeffs0.second()*v,
                coeffs1.first() + coeffs1.second()*v
            ),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalanceModel::eta
(
    const label i,
    const volScalarField::Internal& v
) const
{
    const Pair<dimensionedScalar> coeffs0 = etaCoeffs0(i);
    const Pair<dimensionedScalar> coeffs1 = etaCoeffs1(i);

    return
        max
        (
            min
            (
                coeffs0.first() + coeffs0.second()*v,
                coeffs1.first() + coeffs1.second()*v
            ),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::dimensionedScalar Foam::diameterModels::populationBalanceModel::etaV
(
    const label i,
    const dimensionedScalar& v
) const
{
    const Pair<dimensionedScalar> coeffs0 = etaVCoeffs0(i);
    const Pair<dimensionedScalar> coeffs1 = etaVCoeffs1(i);

    return
        max
        (
            min
            (
                coeffs0.first() + coeffs0.second()/v,
                coeffs1.first() + coeffs1.second()/v
            ),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalanceModel::etaV
(
    const label i,
    const volScalarField::Internal& v
) const
{
    const Pair<dimensionedScalar> coeffs0 = etaVCoeffs0(i);
    const Pair<dimensionedScalar> coeffs1 = etaVCoeffs1(i);

    return
        max
        (
            min
            (
                coeffs0.first() + coeffs0.second()/v,
                coeffs1.first() + coeffs1.second()/v
            ),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::dimensionedScalar Foam::diameterModels::populationBalanceModel::etaV
(
    const labelPair is,
    const distribution& d
) const
{
    // Check that the distribution is generating volumes
    if (d.sampleQ() != 3)
    {
        FatalErrorInFunction
            << "The volumetric allocation coefficient should be evaluated "
            << "with a volumetrically sampled distribution (i.e., sampleQ "
            << "should equal 3)"
            << exit(FatalError);
    }

    // Get the four diameters that bound the sections of the size-group range's
    // basis function. Diameters #1 and #2 are the representative diameters of
    // the first and last size groups in the range. Diameters #0 and #3 are the
    // diameters of the adjacent size groups, or if at an end (or ends), the
    // upper or lower bounding diameter (or diameters) of the distribution.
    //
    //          ^                                                    .
    //          |    o - o - - - - o                                 .
    //  F_first |    .   .         .\                                .
    //          |    .   .         . \                               .
    //          +----+---+---------+--o------------------------------+---->
    //               0   1         2  3                              .
    //          ^    .                                               .
    //          |    .                o - - - - o                    .
    // F_middle |    .               /.         .\                   .
    //          |    .              / .         . \                  .
    //          +----+-------------o--+---------+--o-----------------+---->
    //               .             0  1         2  3                 .
    //          ^    .                                               .
    //          |    .                             o - - - - o - - - o
    //   F_last |    .                            /.         .       .
    //          |    .                           / .         .       .
    //          +----+--------------------------o--+---------+-------+---->
    //               .                          0  1         2       3
    //               .                                               .
    //              dMin                                            dMax
    //
    const bool isFirst = is.first() == 0;
    const bool isLast = is.second() == sizeGroups().size() - 1;
    const scalarField dSphs
    (
        scalarList
        ({
            isFirst
          ? d.min()*(1 - small)
          : sizeGroups()[is.first() - 1].dSph().value(),
            sizeGroups()[is.first()].dSph().value(),
            sizeGroups()[is.second()].dSph().value(),
            isLast
          ? d.max()*(1 + small)
          : sizeGroups()[is.second() + 1].dSph().value()
        })
    );

    // Integrate the distribution, and the distribution divided by volume,
    // across the size-group range's basis function
    const scalarField integralPDFs
    (
        d.integralPDFxPow(dSphs, 0, true)
    );
    const scalarField integralPDFByVs
    (
        d.integralPDFxPow(dSphs, -3, true)
       *6/constant::mathematical::pi
    );

    // Compute the integral of the size-group range's basis function multiplied
    // by the PDF.
    //
    // Between diameters #1 and #2 the basis function is a constant value of
    // one, so the integral is just the same as that of the distribution.
    //
    // Between diameters #0 and #1, and between #2 and #3, the basis function
    // is given by 'C0 + C1/v', where C0 and C1 are thecoefficients given by
    // etaVCoeffs0 and etaVCoeffs1. So, the integral is 'Integral(C0*PDF +
    // C1/v*PDF)'. C0 and C1 are constants, so this can be rearranged to
    // 'C0*Integral(PDF) + C1*Integral(PDF/v)'. The integrals in this
    // expression are those calculated above.
    const Pair<dimensionedScalar> etaVCoeffs0 = this->etaVCoeffs0(is.first());
    const Pair<dimensionedScalar> etaVCoeffs1 = this->etaVCoeffs1(is.second());
    return
        etaVCoeffs0.first().value()
       *(integralPDFs[1] - integralPDFs[0])
      + etaVCoeffs0.second().value()
       *(integralPDFByVs[1] - integralPDFByVs[0])
      + integralPDFs[2] - integralPDFs[1]
      + etaVCoeffs1.first().value()
       *(integralPDFs[3] - integralPDFs[2])
      + etaVCoeffs1.second().value()
       *(integralPDFByVs[3] - integralPDFByVs[2]);
}


Foam::dimensionedScalar Foam::diameterModels::populationBalanceModel::etaV
(
    const label i,
    const distribution& d
) const
{
    const sizeGroup& fi = sizeGroups()[i];

    const PtrList<sizeGroup>& vgSizeGroups = fi.group().sizeGroups();

    // Compute the integral for the size-group's basis function and for the
    // velocity-group's basis function. The allocation coefficient for this
    // size-group is then the ratio of these two integrals. The size-group
    // integral is "scaled up" to be a proportion for the phase, rather than
    // for the population balance as a whole.

    const dimensionedScalar sgIntegralPDFetaV =
        etaV(labelPair(i, i), d);
    const dimensionedScalar vgIntegralPDFetaV =
        etaV(labelPair(vgSizeGroups.first().i(), vgSizeGroups.last().i()), d);

    const dimensionedScalar sgSmall(dimless, rootVSmall/vgSizeGroups.size());
    const dimensionedScalar vgSmall(dimless, rootVSmall);

    return max(sgIntegralPDFetaV, sgSmall)/max(vgIntegralPDFetaV, vgSmall);
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceModel::sigmaWithContinuousPhase
(
    const phaseModel& dispersedPhase
) const
{
    return phaseInterface(dispersedPhase, continuousPhase_).sigma();
}


const Foam::phaseCompressible::momentumTransportModel&
Foam::diameterModels::populationBalanceModel::continuousTurbulence() const
{
    return
        mesh_.lookupType<phaseCompressible::momentumTransportModel>
        (
            continuousPhase_.name()
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalanceModel::Sp(const label i) const
{
    return Sp_[i];
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalanceModel::expansionSu
(
    const label i,
    const UPtrList<const volScalarField>& flds
) const
{
    Pair<tmp<volScalarField::Internal>> tSus = expansionSus(i, flds);

    return
        !tSus.first().valid() ? tSus.second()
      : !tSus.second().valid() ? tSus.first()
      : tSus.first() + tSus.second();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalanceModel::expansionSp(const label i) const
{
    const sizeGroup& fi = sizeGroups()[i];
    const phaseModel& phase = fi.phase();

    tmp<volScalarField::Internal> tSp;

    if (i == 0)
    {
        tSp = negPart(expansionRates_[phase.index()]);
    }
    else
    {
        const sizeGroup& fiMinus1 = sizeGroups()[i - 1];

        tSp =
            negPart(expansionRates_[phase.index()])
           *fi.x()/(fi.x() - fiMinus1.x());
    }

    if (i != sizeGroups().size() - 1)
    {
        const sizeGroup& fiPlus1 = sizeGroups()[i + 1];

        tSp.ref() -=
            posPart(expansionRates_[phase.index()])
           *fi.x()/(fiPlus1.x() - fi.x());
    }
    else
    {
        tSp.ref() += posPart(expansionRates_[phase.index()]);
    }

    return tSp;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalanceModel::modelSourceSu
(
    const label i,
    const UPtrList<const volScalarField>& flds
) const
{
    Pair<tmp<volScalarField::Internal>> tRhoSus = modelSourceRhoSus(i, flds);

    const volScalarField::Internal rho = sizeGroups()[i].phase().rho();

    const dimensionedScalar zeroSu
    (
        (flds.empty() ? dimless : flds[i].dimensions())/dimTime,
        scalar(0)
    );

    return
        tRhoSus.first().valid() && tRhoSus.second().valid()
      ? (tRhoSus.first() + tRhoSus.second())/rho
      : tRhoSus.first().valid() ? tRhoSus.first()/rho
      : tRhoSus.second().valid() ? tRhoSus.second()/rho
      : volScalarField::Internal::New(zeroSu.name(), mesh(), zeroSu);
}


void Foam::diameterModels::populationBalanceModel::solve()
{
    if (!solveOnFinalIterOnly() || fluid_.pimple().finalIter())
    {
        const label nCorr =
            solverDict().lookupBackwardsCompatible<label>
            (
                {"nCorrectors", "nCorr"}
            );

        const scalar tolerance = solverDict().lookup<scalar>("tolerance");

        const bool updateSrc = updateSources();

        if (nCorr > 0 && updateSrc)
        {
            precomputeCoalescenceAndBreakup();
        }
        precomputeExpansion();
        precomputeModelSources();

        int iCorr = 0;
        scalar maxInitialResidual = 1;
        while (++iCorr <= nCorr && maxInitialResidual > tolerance)
        {
            Info<< "populationBalance " << this->name()
                << ": Iteration " << iCorr << endl;

            if (updateSrc)
            {
                computeCoalescenceAndBreakup();
            }
            computeExpansion();
            computeModelSources();

            computeDilatationErrors();

            maxInitialResidual = 0;

            forAll(sizeGroups(), i)
            {
                sizeGroup& fi = sizeGroups_[i];
                const phaseModel& phase = fi.phase();
                const volScalarField& alpha = phase;
                const volScalarField& rho = phase.rho();

                fvScalarMatrix fiEqn
                (
                    fvm::ddt(alpha, fi)
                  + fvm::div(phase.alphaPhi(), fi)
                  + fvm::Sp(-(1 - small)*dilatationErrors_[phase.index()], fi)
                  + fvm::SuSp(-small*dilatationErrors_[phase.index()], fi)
                 ==
                    Su_[i] + fvm::Sp(Sp_[i], fi)
                  + expansionSu(i) + fvm::Sp(expansionSp(i), fi)
                  + modelSourceSu(i)
                  + fluid_.fvModels().source(alpha, rho, fi)/rho
                  - correction
                    (
                        fvm::Sp
                        (
                            max(phase.residualAlpha() - alpha, scalar(0))
                           /mesh().time().deltaT(),
                            fi
                        )
                    )
                );

                fiEqn.relax();

                fluid_.fvConstraints().constrain(fiEqn);

                maxInitialResidual = max
                (
                    fiEqn.solve().initialResidual(),
                    maxInitialResidual
                );

                fluid_.fvConstraints().constrain(fi);
            }
        }

        const volScalarField alphaF0
        (
            sizeGroups().first().phase()*sizeGroups().first()
        );

        const volScalarField alphaFNm1
        (
            sizeGroups().last().phase()*sizeGroups().last()
        );

        Info<< this->name() << " sizeGroup phase fraction first, last = "
            << alphaF0.weightedAverage(mesh().V()).value()
            << ' ' << alphaFNm1.weightedAverage(mesh().V()).value()
            << endl;
    }
}


void Foam::diameterModels::populationBalanceModel::correct()
{
    if (velocityGroupPtrs_.size() <= 1) return;

    // Calculate the total void fraction
    if (alphas_.empty())
    {
        alphas_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", this->name()),
                    fluid_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }
    else
    {
        alphas_() = Zero;
    }

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const phaseModel& phase = iter()->phase();

        alphas_() += max(phase, phase.residualAlpha());
    }

    // Calculate the Sauter mean diameter
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

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const phaseModel& phase = iter()->phase();

        invDsm += max(phase, phase.residualAlpha())/(phase.d()*alphas_());
    }

    if (dsm_.empty())
    {
        dsm_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("d", this->name()),
                    fluid_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                1/tInvDsm
            )
        );
    }
    else
    {
        dsm_() = 1/tInvDsm;
    }

    // Calculate the average velocity
    if (U_.empty())
    {
        U_.set
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("U", this->name()),
                    fluid_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimVelocity, Zero)
            )
        );
    }
    else
    {
        U_() = Zero;
    }

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const phaseModel& phase = iter()->phase();

        U_() += max(phase, phase.residualAlpha())*phase.U()/alphas_();
    }
}


// ************************************************************************* //
