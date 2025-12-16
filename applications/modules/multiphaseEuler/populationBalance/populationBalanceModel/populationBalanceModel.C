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

#include "populationBalance.H"
#include "populationBalanceModel.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "shapeModel.H"
#include "coalescenceModel.H"
#include "daughterSizeDistribution.H"
#include "binary.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "distribution.H"
#include "fvSpecificSource.H"
#include "growthFvScalarFieldSource.H"
#include "oneDimensionalDiscretisation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(populationBalanceModel, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::populationBalanceModel::groupFieldIo
(
    const word& name,
    const label i,
    const phaseModel& phase,
    const IOobject::readOption r,
    const bool registerObject
)
{
    return
        IOobject
        (
            IOobject::groupName
            (
                name + (i == -1 ? "Default" : Foam::name(i)),
                phase.name()
            ),
            phase.mesh().time().name(),
            phase.mesh(),
            r,
            IOobject::AUTO_WRITE,
            registerObject
        );
}


Foam::tmp<Foam::volScalarField> Foam::populationBalanceModel::groupField
(
    const word& name,
    const label i,
    const phaseModel& phase
)
{
    typeIOobject<volScalarField> io
    (
        groupFieldIo(name, i, phase, IOobject::MUST_READ, false)
    );

    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                io.headerOk()
              ? io
              : groupFieldIo(name, -1, phase, IOobject::MUST_READ, false),
                phase.mesh()
            )
        );
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

const Foam::dictionary& Foam::populationBalanceModel::coeffDict() const
{
    return fluid_.optionalSubDict("populationBalanceCoeffs").subDict(name_);
}


void Foam::populationBalanceModel::precomputeCoalescenceAndBreakup()
{
    coalescenceModel_->precompute();

    breakupModel_->precompute();
}


void Foam::populationBalanceModel::birthByCoalescence
(
    const label j,
    const label k,
    const volScalarField::Internal& rate
)
{
    const dimensionedScalar vjk = vs_[j] + vs_[k];

    const volScalarField::Internal alphaFjk
    (
        phases_[j]()*fs_[j]()*phases_[k]()*fs_[k]()
    );

    for (label i = j; i < nGroups(); i++)
    {
        const dimensionedScalar Eta = eta(i, vjk);

        if (Eta.value() == 0) continue;

        volScalarField::Internal Sui
        (
            (j == k ? 0.5 : 1)
           *vs_[i]/(vs_[j]*vs_[k])*Eta*rate*alphaFjk
        );

        Su_[i] += Sui;

        const phaseInterface interfaceij(phases_[i], phases_[j]);

        if (dmdtfs_.found(interfaceij))
        {
            *dmdtfs_[interfaceij] +=
                (interfaceij.index(phases_[i]) == 0 ? +1 : -1)
               *vs_[j]/vjk*Sui*phases_[j].rho();
        }

        const phaseInterface interfaceik(phases_[i], phases_[k]);

        if (dmdtfs_.found(interfaceik))
        {
            *dmdtfs_[interfaceik] +=
                (interfaceik.index(phases_[i]) == 0 ? +1 : -1)
               *vs_[k]/vjk*Sui*phases_[k].rho();
        }

        shapeModel_->addCoalescence(Sui, i, j, k);
    }
}


void Foam::populationBalanceModel::deathByCoalescence
(
    const label i,
    const label j,
    const volScalarField::Internal& rate
)
{
    Sp_[i] -= rate*phases_[i]*fs_[j]*phases_[j]/vs_[j];

    if (i == j) return;

    Sp_[j] -= rate*phases_[j]*phases_[i]*fs_[i]/vs_[i];
}


void Foam::populationBalanceModel::birthByDaughterSizeDistributionBreakup
(
    const label k,
    const volScalarField::Internal& rate
)
{
    for (label i = 0; i <= k; i++)
    {
        const volScalarField::Internal Sui
        (
            rate*phases_[k]*fs_[k]
           *vs_[i]/vs_[k]
           *daughterSizeDistributionBreakupModel_->dsd().nik(i, k)
        );

        Su_[i] += Sui;

        const phaseInterface interface(phases_[i], phases_[k]);

        if (dmdtfs_.found(interface))
        {
            *dmdtfs_[interface] +=
                (interface.index(phases_[i]) == 0 ? +1 : -1)
               *Sui*phases_[k].rho();
        }

        shapeModel_->addBreakup(Sui, i, k);
    }
}


void Foam::populationBalanceModel::deathByDaughterSizeDistributionBreakup
(
    const label i,
    const volScalarField::Internal& rate
)
{
    Sp_[i] -= rate*phases_[i];
}


void Foam::populationBalanceModel::birthByBinaryBreakup
(
    const label i,
    const label j,
    const volScalarField::Internal& rate
)
{
    const volScalarField::Internal Su(rate*phases_[j]*fs_[j]);

    {
        const volScalarField::Internal Sui
        (
            vs_[i]*binaryBreakupDeltas_[i][j]/vs_[j]*Su
        );

        Su_[i] += Sui;

        const phaseInterface interfaceij(phases_[i], phases_[j]);

        if (dmdtfs_.found(interfaceij))
        {
            *dmdtfs_[interfaceij] +=
                (interfaceij.index(phases_[i]) == 0 ? +1 : -1)
               *Sui*phases_[j].rho();
        }

        shapeModel_->addBreakup(Sui, i, j);
    }

    const dimensionedScalar v = vs_[j] - vs_[i];

    for (label k = 0; k <= j; k ++)
    {
        const dimensionedScalar Eta = eta(k, v);

        if (Eta.value() == 0) continue;

        const volScalarField::Internal Suk
        (
            vs_[k]*binaryBreakupDeltas_[i][j]*Eta/vs_[j]*Su
        );

        Su_[k] += Suk;

        const phaseInterface interfacekj(phases_[k], phases_[j]);

        if (dmdtfs_.found(interfacekj))
        {
            *dmdtfs_[interfacekj] +=
                (interfacekj.index(phases_[k]) == 0 ? +1 : -1)
               *Suk*phases_[j].rho();
        }

        shapeModel_->addBreakup(Suk, k, j);
    }
}


void Foam::populationBalanceModel::deathByBinaryBreakup
(
    const label j,
    const label i,
    const volScalarField::Internal& rate
)
{
    Sp_[i] -= rate*phases_[i]*binaryBreakupDeltas_[j][i];
}


void Foam::populationBalanceModel::computeCoalescenceAndBreakup()
{
    shapeModel_->reset();

    forAll(fs_, i)
    {
        Su_[i] = Zero;
        Sp_[i] = Zero;
    }

    forAllIter(dmdtfTable, dmdtfs_, dmdtfIter)
    {
        *dmdtfIter() = Zero;
    }

    if (coalescenceModel_->coalesces())
    {
        forAll(coalescencePairs_, coalescencePairi)
        {
            const label i = coalescencePairs_[coalescencePairi].first();
            const label j = coalescencePairs_[coalescencePairi].second();

            tmp<volScalarField::Internal> trate = coalescenceModel_->rate(i, j);

            birthByCoalescence(i, j, trate());

            deathByCoalescence(i, j, trate());
        }
    }

    if (daughterSizeDistributionBreakupModel_)
    {
        forAll(fs_, i)
        {
            tmp<volScalarField::Internal> trate =
                daughterSizeDistributionBreakupModel_->rate(i);

            birthByDaughterSizeDistributionBreakup(i, trate());

            deathByDaughterSizeDistributionBreakup(i, trate());
        }
    }

    if (binaryBreakupModel_)
    {
        forAll(binaryBreakupPairs_, binaryBreakupPairi)
        {
            const label i = binaryBreakupPairs_[binaryBreakupPairi].first();
            const label j = binaryBreakupPairs_[binaryBreakupPairi].second();

            tmp<volScalarField::Internal> trate =
                binaryBreakupModel_->rate(j, i);

            birthByBinaryBreakup(j, i, trate());

            deathByBinaryBreakup(j, i, trate());
        }
    }
}


void Foam::populationBalanceModel::precomputeExpansion()
{
    forAll(uniquePhases_, uniquePhasei)
    {
        const phaseModel& phase = uniquePhases_[uniquePhasei];
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
Foam::populationBalanceModel::expansionSus
(
    const label i,
    const UPtrList<volScalarField>& flds
) const
{
    auto fiFld = [&](const label deltai)
    {
        return
            flds.empty()
          ? tmp<volScalarField::Internal>(fs_[i + deltai])
          : fs_[i + deltai]()*flds[i + deltai]();
    };

    Pair<tmp<volScalarField::Internal>> tSus;

    if (i != 0)
    {
        tSus.first() =
            posPart(expansionRates_[phases_[i - 1].index()])
           *vs_[i]/(vs_[i] - vs_[i - 1])
           *fiFld(-1);
    }

    if (i != nGroups() - 1)
    {
        tSus.second() =
          - negPart(expansionRates_[phases_[i + 1].index()])
           *vs_[i]/(vs_[i + 1] - vs_[i])
           *fiFld(+1);
    }

    return tSus;
}


void Foam::populationBalanceModel::computeExpansion()
{
    forAllIter(dmdtfTable, expansionDmdtfs_, dmdtfIter)
    {
        *dmdtfIter() = Zero;
    }

    for
    (
        label uniquePhasei = 0;
        uniquePhasei < uniquePhases_.size() - 1;
        ++ uniquePhasei
    )
    {
        const phaseModel& phase0 = uniquePhases_[uniquePhasei];
        const phaseModel& phase1 = uniquePhases_[uniquePhasei + 1];

        const label i0 = uniqueDiameters_[uniquePhasei].iLast();
        const label i1 = uniqueDiameters_[uniquePhasei + 1].iFirst();

        Pair<tmp<volScalarField::Internal>> tSus0 = expansionSus(i0);
        Pair<tmp<volScalarField::Internal>> tSus1 = expansionSus(i1);

        const phaseInterface interface01(phase0, phase1);

        *expansionDmdtfs_[interface01] +=
            (interface01.index(phase0) == 0 ? -1 : +1)
           *(- tSus0.second()*phase0.rho() + tSus1.first()*phase1.rho());
    }
}


void Foam::populationBalanceModel::precomputeModelSources()
{
    // nothing to do
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::populationBalanceModel::modelSourceRhoSus
(
    const label i,
    const UPtrList<volScalarField>& flds
) const
{
    const volScalarField& fldi = flds.empty() ? fs_[i] : flds[i];

    Pair<tmp<volScalarField::Internal>> tRhoSus;

    forAll(tRhoSus, Sui)
    {
        const label iPopBalEnd = Sui == 0 ? 0 : nGroups() - 1;

        if (i == iPopBalEnd) continue;

        const label iVelGrpEnd =
            Sui == 0 ? diameters_[i].iFirst() : diameters_[i].iLast();

        if (i != iVelGrpEnd) continue;

        const label iOther = i + (Sui == 0 ? -1 : +1);

        forAll(fluid_.fvModels(), modeli)
        {
            if (!isA<fvSpecificSource>(fluid_.fvModels()[modeli])) continue;

            const fvSpecificSource& source =
                refCast<const fvSpecificSource>(fluid_.fvModels()[modeli]);

            if
            (
                source.addsSupToField(phases_[iOther].volScalarField::name())
             && isA<growthFvScalarFieldSource>(fldi.sources()[source.name()])
            )
            {
                const growthFvScalarFieldSource& growthSource =
                    refCast<const growthFvScalarFieldSource>
                    (
                        fldi.sources()[source.name()]
                    );

                const volScalarField::Internal S(source.S(fs_[iOther].name()));

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


void Foam::populationBalanceModel::computeModelSources()
{
    forAllIter(dmdtfTable, modelSourceDmdtfs_, dmdtfIter)
    {
        *dmdtfIter() = Zero;
    }

    forAll(uniquePhases_, uniquePhasei)
    {
        const label i0 = uniqueDiameters_[uniquePhasei].iLast();

        if (i0 == nGroups() - 1) continue;

        const label i1 = i0 + 1;

        Pair<tmp<volScalarField::Internal>> tRhoSus0 = modelSourceRhoSus(i0);
        Pair<tmp<volScalarField::Internal>> tRhoSus1 = modelSourceRhoSus(i1);

        const phaseInterface interface01(phases_[i0], phases_[i1]);
        const scalar sign = interface01.index(phases_[i0]) == 0 ? -1 : +1;

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


void Foam::populationBalanceModel::computeDilatationErrors()
{
    PtrList<volScalarField::Internal> modelSourceDmdts(fluid_.phases().size());
    forAllConstIter(dmdtfTable, modelSourceDmdtfs_, dmdtfIter)
    {
        const phaseInterface interface(fluid_, dmdtfIter.key());

        addField(interface.phase1(), "dmdt", *dmdtfIter(), modelSourceDmdts);
        addField(interface.phase2(), "dmdt", - *dmdtfIter(), modelSourceDmdts);
    }

    forAll(uniquePhases_, uniquePhasei)
    {
        const phaseModel& phase = uniquePhases_[uniquePhasei];
        const diameterModels::populationBalance& diameter =
            uniqueDiameters_[uniquePhasei];
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        dilatationErrors_.set
        (
            phase.index(),
            fvc::ddt(alpha)()() + fvc::div(phase.alphaPhi())()()
          - (fluid_.fvModels().source(alpha, rho) & rho)()()/rho()
        );

        for (label i = diameter.iFirst(); i <= diameter.iLast(); ++ i)
        {
            dilatationErrors_[phase.index()] -=
                Su_[i] + expansionSu(i) + (Sp_[i] + expansionSp(i))*fs_[i];
        }

        if (modelSourceDmdts.set(phase.index()))
        {
            dilatationErrors_[phase.index()] -=
                modelSourceDmdts[phase.index()]/rho;
        }
    }
}


bool Foam::populationBalanceModel::updateSources()
{
    const bool result = sourceUpdateCounter_ % sourceUpdateInterval() == 0;

    ++ sourceUpdateCounter_;

    return result;
}


Foam::Pair<Foam::dimensionedScalar> Foam::populationBalanceModel::etaCoeffs0
(
    const label i
) const
{
    return
        i == 0
      ? Pair<dimensionedScalar>
        (
            dimensionedScalar(dimless, scalar(0)),
            1/vs_[i]
        )
      : Pair<dimensionedScalar>
        (
            -vs_[i-1]/(vs_[i] - vs_[i-1]),
            1/(vs_[i] - vs_[i-1])
        );
}


Foam::Pair<Foam::dimensionedScalar> Foam::populationBalanceModel::etaCoeffs1
(
    const label i
) const
{
    return
        i == vs_.size() - 1
      ? Pair<dimensionedScalar>
        (
            dimensionedScalar(dimless, scalar(0)),
            1/vs_[i]
        )
      : Pair<dimensionedScalar>
        (
            vs_[i+1]/(vs_[i+1] - vs_[i]),
            -1/(vs_[i+1] - vs_[i])
        );
}


Foam::Pair<Foam::dimensionedScalar> Foam::populationBalanceModel::etaVCoeffs0
(
    const label i
) const
{
    return
        i == 0
      ? Pair<dimensionedScalar>
        (
            dimensionedScalar(dimless, scalar(1)),
            dimensionedScalar(dimVolume, scalar(0))
        )
      : Pair<dimensionedScalar>
        (
            -1/vs_[i-1]/(1/vs_[i] - 1/vs_[i-1]),
            1/(1/vs_[i] - 1/vs_[i-1])
        );
}


Foam::Pair<Foam::dimensionedScalar> Foam::populationBalanceModel::etaVCoeffs1
(
    const label i
) const
{
    return
        i == vs_.size() - 1
      ? Pair<dimensionedScalar>
        (
            dimensionedScalar(dimless, scalar(1)),
            dimensionedScalar(dimVolume, scalar(0))
        )
      : Pair<dimensionedScalar>
        (
            1/vs_[i+1]/(1/vs_[i+1] - 1/vs_[i]),
            -1/(1/vs_[i+1] - 1/vs_[i])
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceModel::populationBalanceModel
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
    phases_(),
    uniquePhases_(),
    diameters_(),
    uniqueDiameters_(),
    fs_(),
    dSphs_(),
    vs_(),
    Su_(),
    Sp_(),
    dmdtfs_(),
    expansionDmdtfs_(),
    modelSourceDmdtfs_(),
    expansionRates_(fluid_.phases().size()),
    dilatationErrors_(fluid_.phases().size()),
    shapeModel_(nullptr),
    coalescenceModel_(),
    coalescencePairs_(),
    breakupModel_(),
    daughterSizeDistributionBreakupModel_(nullptr),
    binaryBreakupModel_(nullptr),
    binaryBreakupDeltas_(),
    binaryBreakupPairs_(),
    alphas_(),
    dsm_(),
    U_(),
    sourceUpdateCounter_(0)
{
    Info<< "Population balance model: " << name << incrIndent << endl;

    // Build the phase-reference lists
    for
    (
        label fluidPhasei = 0, nGroups = 0;
        fluidPhasei < fluid.phases().size();
        ++ fluidPhasei
    )
    {
        const phaseModel& phase = fluid.phases()[fluidPhasei];

        if (!isA<diameterModels::populationBalance>(phase.diameter())) continue;

        const diameterModels::populationBalance& diameter =
            refCast<const diameterModels::populationBalance>(phase.diameter());

        if (diameter.popBalName() != name_) continue;

        uniquePhases_.append(&phase);
        uniqueDiameters_.append(&diameter);

        phases_.resize(nGroups + diameter.nGroups());
        diameters_.resize(nGroups + diameter.nGroups());
        for (label i = 0; i < diameter.nGroups(); ++ i)
        {
            phases_.set(nGroups + i, &phase);
            diameters_.set(nGroups + i, &diameter);
        }

        nGroups += diameter.nGroups();
    }

    // Check that there are sufficiently many groups
    if (nGroups() < 3)
    {
        FatalErrorInFunction
            << "A population balance model requires a minimum of 3 groups but "
            << name_ << " only has " << nGroups() << exit(FatalError);
    }

    // Read/initialise the group fraction fields
    fs_.resize(nGroups());
    forAll(fs_, i)
    {
        fs_.set
        (
            i,
            new volScalarField
            (
                groupFieldIo("f", i, phases_[i]),
                groupField("f", i, phases_[i])
            )
        );
    }

    // Create a discretisation in representative spherical diameter space
    dSphs_ =
        oneDimensionalDiscretisation::New
        (
            "dSph",
            dimLength,
            nGroups(),
            coeffDict().subDict("sphericalDiameters")
        )->dimensionedCoordinates();

    // Build the groups' representative volumes
    vs_.resize(nGroups());
    forAll(fs_, i)
    {
        vs_.set
        (
            i,
            new dimensionedScalar
            (
                "v" + Foam::name(i),
                constant::mathematical::pi/6*pow3(dSphs_[i])
            )
        );
    }

    // Print some confirmatory information about the setup
    forAll(uniquePhases_, uniquePhasei)
    {
        const phaseModel& phase = uniquePhases_[uniquePhasei];
        const diameterModels::populationBalance& diameter =
            uniqueDiameters_[uniquePhasei];

        Info<< indent << "Phase: " << phase.name() << incrIndent << endl;

        for (label i = diameter.iFirst(); i <= diameter.iLast(); ++ i)
        {
            Info<< indent << "Group #" << i
                << ": dSph = " << dSphs_[i].value()
                << ", min/average/max fraction = "
                << min(fs_[i]()).value() << '/'
                << average(fs_[i]()).value() << '/'
                << max(fs_[i]()).value() << endl;
        }

        Info<< decrIndent;
    }

    // Create source terms
    Su_.setSize(nGroups());
    Sp_.setSize(nGroups());
    forAll(Su_, i)
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
    forAll(uniquePhases_, uniquePhasei)
    {
        for
        (
            label uniquePhasej = 0;
            uniquePhasej < uniquePhasei;
            ++ uniquePhasej
        )
        {
            const phaseInterface interface
            (
                uniquePhases_[uniquePhasei],
                uniquePhases_[uniquePhasej]
            );

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

    using namespace populationBalance;

    // Select the shape model
    shapeModel_.set(shapeModel::New(coeffDict(), *this).ptr());

    // Select coalescence model
    coalescenceModel_.set(coalescenceModel::New(*this, coeffDict()).ptr());
    if (coalescenceModel_->coalesces())
    {
        forAll(fs_, i)
        {
            for (label j = 0; j <= i; j++)
            {
                coalescencePairs_.append(labelPair(i, j));
            }
        }
    }

    // Select breakup model
    breakupModel_.set(breakupModel::New(*this, coeffDict()).ptr());
    if (isA<breakupModels::daughterSizeDistribution>(breakupModel_()))
    {
        daughterSizeDistributionBreakupModel_ =
            &refCast<breakupModels::daughterSizeDistribution>(breakupModel_());
    }
    if (isA<breakupModels::binary>(breakupModel_()))
    {
        binaryBreakupModel_ =
            &refCast<breakupModels::binary>(breakupModel_());
    }
    if (binaryBreakupModel_)
    {
        binaryBreakupDeltas_.setSize(nGroups());

        forAll(fs_, i)
        {
            binaryBreakupDeltas_.set
            (
                i,
                new PtrList<dimensionedScalar>(nGroups())
            );

            const dimensionedScalar vMid0
            (
                i == 0 ? vs_.first() : (vs_[i - 1] + vs_[i])/2
            );
            const dimensionedScalar vMid1
            (
                i == nGroups() - 1 ? vs_.last() : (vs_[i] + vs_[i + 1])/2
            );

            forAll(fs_, j)
            {
                binaryBreakupDeltas_[i].set
                (
                    j,
                    new dimensionedScalar
                    (
                        "binaryBreakupDelta_"
                      + Foam::name(i)
                      + "_"
                      + Foam::name(j),
                        vMid0 < 0.5*vs_[j] && 0.5*vs_[j] < vMid1
                      ? mag(0.5*vs_[j] - vMid0)
                      : 0.5*vs_[j] < vMid0
                      ? dimensionedScalar(dimVolume, scalar(0))
                      : vMid1 - vMid0
                    )
                );
            }
        }

        forAll(fs_, i)
        {
            label j = 0;

            while (binaryBreakupDeltas_[j][i].value() != 0)
            {
                binaryBreakupPairs_.append(labelPair(i, j));
                j++;
            }
        }
    }

    // Everything is now available so correct the diameter models and the
    // population balance model so that all publicly available data is valid
    forAll(uniquePhases_, uniquePhasei)
    {
        const_cast<diameterModels::populationBalance&>
        (
            uniqueDiameters_[uniquePhasei]
        ).correct();
    }

    correct();

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceModel::~populationBalanceModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::populationBalance::shapeModel&
Foam::populationBalanceModel::shape() const
{
    return shapeModel_();
}


Foam::tmp<Foam::volScalarField> Foam::populationBalanceModel::a
(
    const label i
) const
{
    return shapeModel_->a(i);
}


Foam::tmp<Foam::volScalarField> Foam::populationBalanceModel::d
(
    const label i
) const
{
    return shapeModel_->d(i);
}


Foam::dimensionedScalar Foam::populationBalanceModel::eta
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


Foam::tmp<Foam::volScalarField::Internal> Foam::populationBalanceModel::eta
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


Foam::dimensionedScalar Foam::populationBalanceModel::etaV
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
Foam::populationBalanceModel::etaV
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


Foam::dimensionedScalar Foam::populationBalanceModel::etaV
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

    // Get the four diameters that bound the sections of the group range's
    // basis function. Diameters #1 and #2 are the representative diameters of
    // the first and last groups in the range. Diameters #0 and #3 are the
    // diameters of the adjacent groups, or if at an end (or ends), the upper
    // or lower bounding diameter (or diameters) of the distribution.
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
    const bool isLast = is.second() == nGroups() - 1;
    const scalarField dSphs
    (
        scalarList
        ({
            isFirst
          ? d.min()*(1 - small)
          : dSphs_[is.first() - 1].value(),
            dSphs_[is.first()].value(),
            dSphs_[is.second()].value(),
            isLast
          ? d.max()*(1 + small)
          : dSphs_[is.second() + 1].value()
        })
    );

    // Integrate the distribution, and the distribution divided by volume,
    // across the group range's basis function
    const scalarField integralPDFs
    (
        d.integralPDFxPow(dSphs, 0, true)
    );
    const scalarField integralPDFByVs
    (
        d.integralPDFxPow(dSphs, -3, true)
       *6/constant::mathematical::pi
    );

    // Compute the integral of the group range's basis function multiplied by
    // the PDF.
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


Foam::dimensionedScalar Foam::populationBalanceModel::etaV
(
    const label i,
    const distribution& d
) const
{
    const label i0 = diameters_[i].iFirst();
    const label iNm1 = diameters_[i].iLast();

    // Compute the integral for the group's basis function and for the phase's
    // basis function. The allocation coefficient for this group is then the
    // ratio of these two integrals. This means the group integral is "scaled
    // up" to be a proportion for the phase, rather than for the population
    // balance as a whole.

    const dimensionedScalar groupIntegralPDFetaV = etaV(labelPair(i, i), d);
    const dimensionedScalar phaseIntegralPDFetaV = etaV(labelPair(i0, iNm1), d);

    static const dimensionedScalar dimlessRootVSmall(dimless, rootVSmall);

    return
        max(groupIntegralPDFetaV, dimlessRootVSmall/diameters_[i].nGroups())
       /max(phaseIntegralPDFetaV, dimlessRootVSmall);
}


const Foam::tmp<Foam::volScalarField>
Foam::populationBalanceModel::sigmaWithContinuousPhase
(
    const phaseModel& dispersedPhase
) const
{
    return phaseInterface(dispersedPhase, continuousPhase_).sigma();
}


const Foam::tmp<Foam::volScalarField>
Foam::populationBalanceModel::sigmaWithContinuousPhase
(
    const label i
) const
{
    return phaseInterface(phases_[i], continuousPhase_).sigma();
}


const Foam::phaseCompressible::momentumTransportModel&
Foam::populationBalanceModel::continuousTurbulence() const
{
    return
        mesh_.lookupType<phaseCompressible::momentumTransportModel>
        (
            continuousPhase_.name()
        );
}


Foam::tmp<Foam::volScalarField::Internal> Foam::populationBalanceModel::Sp
(
    const label i
) const
{
    return Sp_[i];
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalanceModel::expansionSu
(
    const label i,
    const UPtrList<volScalarField>& flds
) const
{
    Pair<tmp<volScalarField::Internal>> tSus = expansionSus(i, flds);

    return
        !tSus.first().valid() ? tSus.second()
      : !tSus.second().valid() ? tSus.first()
      : tSus.first() + tSus.second();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalanceModel::expansionSp
(
    const label i
) const
{
    const volScalarField::Internal& expansionRate =
        expansionRates_[phases_[i].index()];

    tmp<volScalarField::Internal> tSp;

    if (i == 0)
    {
        tSp = negPart(expansionRate);
    }
    else
    {
        tSp = negPart(expansionRate)*vs_[i]/(vs_[i] - vs_[i - 1]);
    }

    if (i != nGroups() - 1)
    {
        tSp.ref() -= posPart(expansionRate)*vs_[i]/(vs_[i + 1] - vs_[i]);
    }
    else
    {
        tSp.ref() += posPart(expansionRate);
    }

    return tSp;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalanceModel::modelSourceSu
(
    const label i,
    const UPtrList<volScalarField>& flds
) const
{
    Pair<tmp<volScalarField::Internal>> tRhoSus = modelSourceRhoSus(i, flds);

    const dimensionedScalar zeroSu
    (
        (flds.empty() ? dimless : flds[i].dimensions())/dimTime,
        scalar(0)
    );

    return
        tRhoSus.first().valid() && tRhoSus.second().valid()
      ? (tRhoSus.first() + tRhoSus.second())/phases_[i].rho()
      : tRhoSus.first().valid() ? tRhoSus.first()/phases_[i].rho()
      : tRhoSus.second().valid() ? tRhoSus.second()/phases_[i].rho()
      : volScalarField::Internal::New(zeroSu.name(), mesh(), zeroSu);
}


void Foam::populationBalanceModel::solve()
{
    if (solveOnFinalIterOnly() && !fluid_.pimple().finalIter())
    {
        return;
    }

    const label nCorr =
        solverDict().lookupBackwardsCompatible<label>({"nCorrectors", "nCorr"});

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

        forAll(fs_, i)
        {
            volScalarField& fi = fs_[i];
            const phaseModel& phase = phases_[i];
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

    const volScalarField alphaF0(phases_.first()*fs_.first());
    const volScalarField alphaFNm1(phases_.last()*fs_.last());

    Info<< "populationBalance " << this->name() << ": Group fraction "
        << "first/last = " << alphaF0.weightedAverage(mesh().V()).value()
        << '/' << alphaFNm1.weightedAverage(mesh().V()).value() << endl;

    if (solverDict().lookupOrDefault<Switch>("scale", true))
    {
        Info<< "populationBalance " << this->name()
            << ": Scaling group fractions " << endl;

        forAll(fs_, i)
        {
            fs_[i].max(0);
        }

        forAll(uniquePhases_, uniquePhasei)
        {
            const diameterModels::populationBalance& diameter =
                uniqueDiameters_[uniquePhasei];

            const volScalarField::Internal fSum(diameter.fSum());

            for (label i = diameter.iFirst(); i <= diameter.iLast(); ++ i)
            {
                fs_[i].internalFieldRef() /= fSum;

                fs_[i].correctBoundaryConditions();
            }
        }
    }
    else
    {
        forAll(uniquePhases_, uniquePhasei)
        {
            const diameterModels::populationBalance& diameter =
                uniqueDiameters_[uniquePhasei];

            const volScalarField::Internal fSum(diameter.fSum());

            Info<< diameter.phase().name()
                << ": Group fraction sum min/average/max = "
                << min(fSum).value() << '/'
                << fSum.weightedAverage(mesh().V()).value() << '/'
                << max(fSum).value() << endl;
        }
    }

    shapeModel_->solve();
}


void Foam::populationBalanceModel::correct()
{
    shapeModel_->correct();

    if (uniquePhases_.size() <= 1) return;

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

    forAll(uniquePhases_, uniquePhasei)
    {
        alphas_() +=
            max
            (
                uniquePhases_[uniquePhasei],
                uniquePhases_[uniquePhasei].residualAlpha()
            );
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

    forAll(uniquePhases_, uniquePhasei)
    {
        invDsm +=
            max
            (
                uniquePhases_[uniquePhasei],
                uniquePhases_[uniquePhasei].residualAlpha()
            )
           /alphas_()
           /uniquePhases_[uniquePhasei].d();
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

    forAll(uniquePhases_, uniquePhasei)
    {
        U_() +=
            max
            (
                uniquePhases_[uniquePhasei],
                uniquePhases_[uniquePhasei].residualAlpha()
            )
           /alphas_()
           *uniquePhases_[uniquePhasei].U();
    }
}


bool Foam::populationBalanceModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
