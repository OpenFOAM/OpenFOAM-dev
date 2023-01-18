/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
#include "coalescenceModel.H"
#include "breakupModel.H"
#include "binaryBreakupModel.H"
#include "driftModel.H"
#include "nucleationModel.H"
#include "interfaceSurfaceTensionModel.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSup.H"
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
                velocityGroupPtrs_.insert(velGroup.phase().name(), &velGroup);

                dilatationErrors_.insert
                (
                    velGroup.phase().name(),
                    volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "dilatationError",
                                velGroup.phase().name()
                            ),
                            fluid_.time().name(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedScalar(inv(dimTime), 0)
                    )
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
        sizeGroups().size() != 0
        &&
        group.x().value() <= sizeGroups().last().x().value()
    )
    {
        FatalErrorInFunction
            << "Size groups must be entered according to their representative"
            << " size"
            << exit(FatalError);
    }

    sizeGroups_.resize(sizeGroups().size() + 1);
    sizeGroups_.set(sizeGroups().size() - 1, &group);

    if (sizeGroups().size() == 1)
    {
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups().last().x()
            )
        );

        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups().last().x()
            )
        );
    }
    else
    {
        v_.last() =
            0.5
           *(
                sizeGroups()[sizeGroups().size()-2].x()
              + sizeGroups().last().x()
            );

        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups().last().x()
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
                fluid_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(inv(dimTime), 0)
        )
    );

    Sp_.append
    (
        new volScalarField
        (
            IOobject
            (
                "Sp",
                fluid_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(inv(dimTime), 0)
        )
    );
}


void Foam::diameterModels::populationBalanceModel::initialiseDmdtfs()
{
    forAllConstIter
    (
        HashTable<const diameterModels::velocityGroup*>,
        velocityGroupPtrs(),
        iter1
    )
    {
        const diameterModels::velocityGroup& velGrp1 = *iter1();

        forAllConstIter
        (
            HashTable<const diameterModels::velocityGroup*>,
            velocityGroupPtrs(),
            iter2
        )
        {
            const diameterModels::velocityGroup& velGrp2 = *iter2();

            const phaseInterface interface(velGrp1.phase(), velGrp2.phase());

            if
            (
                &velGrp1 != &velGrp2
                &&
                !dmdtfs_.found(interface)
            )
            {
                fluid_.validateMassTransfer
                    <diameterModels::populationBalanceModel>(interface);

                dmdtfs_.insert
                (
                    interface,
                    new volScalarField
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
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::precompute()
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

    forAll(drift_, model)
    {
        drift_[model].precompute();
    }

    forAll(nucleation_, model)
    {
        nucleation_[model].precompute();
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

    dimensionedScalar Eta;
    dimensionedScalar v = fj.x() + fk.x();

    for (label i = j; i < sizeGroups().size(); i++)
    {
        Eta = eta(i, v);

        if (Eta.value() == 0) continue;

        const sizeGroup& fi = sizeGroups()[i];

        if (j == k)
        {
            Sui_ =
                0.5*fi.x()/(fj.x()*fk.x())*Eta
               *coalescenceRate_()*fj*fj.phase()*fk*fk.phase();
        }
        else
        {
            Sui_ =
                fi.x()/(fj.x()*fk.x())*Eta
               *coalescenceRate_()*fj*fj.phase()*fk*fk.phase();
        }

        Su_[i] += Sui_;

        const phaseInterface interfaceij(fi.phase(), fj.phase());

        if (dmdtfs_.found(interfaceij))
        {
            const scalar dmdtSign =
                interfaceij.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfaceij] += dmdtSign*fj.x()/v*Sui_*fj.phase().rho();
        }

        const phaseInterface interfaceik(fi.phase(), fk.phase());

        if (dmdtfs_.found(interfaceik))
        {
            const scalar dmdtSign =
                interfaceik.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfaceik] += dmdtSign*fk.x()/v*Sui_*fk.phase().rho();
        }

        sizeGroups_[i].shapeModelPtr()->addCoalescence(Sui_, fj, fk);
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

    Sp_[i] += coalescenceRate_()*fi.phase()*fj*fj.phase()/fj.x();

    if (i != j)
    {
        Sp_[j] += coalescenceRate_()*fj.phase()*fi*fi.phase()/fi.x();
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

        Sui_ =
            fi.x()*breakupModels_[model].dsdPtr()().nik(i, k)/fk.x()
           *breakupRate_()*fk*fk.phase();

        Su_[i] += Sui_;

        const phaseInterface interface(fi.phase(), fk.phase());

        if (dmdtfs_.found(interface))
        {
            const scalar dmdtSign =
                interface.index(fi.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interface] += dmdtSign*Sui_*fk.phase().rho();
        }

        sizeGroups_[i].shapeModelPtr()->addBreakup(Sui_, fk);
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBreakup(const label i)
{
    Sp_[i] += breakupRate_()*sizeGroups()[i].phase();
}


void Foam::diameterModels::populationBalanceModel::calcDeltas()
{
    forAll(sizeGroups(), i)
    {
        if (delta_[i].empty())
        {
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
                 &&
                    0.5*fj.x().value() < v_[i+1].value()
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
}


void Foam::diameterModels::populationBalanceModel::birthByBinaryBreakup
(
    const label i,
    const label j
)
{
    const sizeGroup& fi = sizeGroups()[i];
    const sizeGroup& fj = sizeGroups()[j];

    const volScalarField Su(binaryBreakupRate_()*fj*fj.phase());

    Sui_ = fi.x()*delta_[i][j]/fj.x()*Su;

    Su_[i] += Sui_;

    sizeGroups_[i].shapeModelPtr()->addBreakup(Sui_, fj);

    const phaseInterface interfaceij(fi.phase(), fj.phase());

    if (dmdtfs_.found(interfaceij))
    {
        const scalar dmdtSign =
            interfaceij.index(fi.phase()) == 0 ? +1 : -1;

        *dmdtfs_[interfaceij] += dmdtSign*Sui_*fj.phase().rho();
    }

    dimensionedScalar Eta;
    dimensionedScalar v = fj.x() - fi.x();

    for (label k = 0; k <= j; k++)
    {
        Eta = eta(k, v);

        if (Eta.value() == 0) continue;

        const sizeGroup& fk = sizeGroups()[k];

        volScalarField& Suk = Sui_;

        Suk = fk.x()*delta_[i][j]*Eta/fj.x()*Su;

        Su_[k] += Suk;

        const phaseInterface interfacekj(fk.phase(), fj.phase());

        if (dmdtfs_.found(interfacekj))
        {
            const scalar dmdtSign =
                interfacekj.index(fk.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfacekj] += dmdtSign*Suk*fj.phase().rho();
        }

        sizeGroups_[k].shapeModelPtr()->addBreakup(Suk, fj);
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBinaryBreakup
(
    const label j,
    const label i
)
{
    Sp_[i] += sizeGroups()[i].phase()*binaryBreakupRate_()*delta_[j][i];
}


void Foam::diameterModels::populationBalanceModel::drift
(
    const label i,
    driftModel& model
)
{
    model.addToDriftRate(driftRate_(), i);

    const sizeGroup& fp = sizeGroups()[i];

    if (i < sizeGroups().size() - 1)
    {
        const sizeGroup& fe = sizeGroups()[i+1];
        volScalarField& Sue = Sui_;

        Sp_[i] += 1/(fe.x() - fp.x())*pos(driftRate_())*driftRate_();

        Sue =
            fe.x()/(fp.x()*(fe.x() - fp.x()))*pos(driftRate_())*driftRate_()*fp;

        Su_[i+1] += Sue;

        const phaseInterface interfacepe(fp.phase(), fe.phase());

        if (dmdtfs_.found(interfacepe))
        {
            const scalar dmdtSign =
                interfacepe.index(fp.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfacepe] -= dmdtSign*Sue*fp.phase().rho();
        }

        sizeGroups_[i+1].shapeModelPtr()->addDrift(Sue, fp, model);
    }

    if (i == sizeGroups().size() - 1)
    {
        Sp_[i] -= pos(driftRate_())*driftRate_()/fp.x();
    }

    if (i > 0)
    {
        const sizeGroup& fw = sizeGroups()[i-1];
        volScalarField& Suw = Sui_;

        Sp_[i] += 1/(fw.x() - fp.x())*neg(driftRate_())*driftRate_();

        Suw =
            fw.x()/(fp.x()*(fw.x() - fp.x()))*neg(driftRate_())*driftRate_()*fp;

        Su_[i-1] += Suw;

        const phaseInterface interfacepw(fp.phase(), fw.phase());

        if (dmdtfs_.found(interfacepw))
        {
            const scalar dmdtSign =
                interfacepw.index(fp.phase()) == 0 ? +1 : -1;

            *dmdtfs_[interfacepw] -= dmdtSign*Suw*fp.phase().rho();
        }

        sizeGroups_[i-1].shapeModelPtr()->addDrift(Suw, fp, model);
    }

    if (i == 0)
    {
        Sp_[i] -= neg(driftRate_())*driftRate_()/fp.x();
    }
}


void Foam::diameterModels::populationBalanceModel::nucleation
(
    const label i,
    nucleationModel& model
)
{
    const sizeGroup& fi = sizeGroups()[i];

    model.addToNucleationRate(nucleationRate_(), i);

    Sui_ = fi.x()*nucleationRate_();

    Su_[i] += Sui_;

    sizeGroups_[i].shapeModelPtr()->addNucleation(Sui_, fi, model);
}


void Foam::diameterModels::populationBalanceModel::sources()
{
    forAll(sizeGroups(), i)
    {
        sizeGroups_[i].shapeModelPtr()->reset();
        Su_[i] = Zero;
        Sp_[i] = Zero;
    }

    forAllIter(phaseSystem::dmdtfTable, dmdtfs_, pDmdtIter)
    {
        *pDmdtIter() = Zero;
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

    forAll(sizeGroups(), i)
    {
        forAll(breakupModels_, model)
        {
            breakupModels_[model].setBreakupRate(breakupRate_(), i);

            birthByBreakup(i, model);

            deathByBreakup(i);
        }

        forAll(drift_, model)
        {
            driftRate_() = Zero;

            drift(i, drift_[model]);
        }

        forAll(nucleation_, model)
        {
            nucleationRate_() = Zero;

            nucleation(i, nucleation_[model]);
        }
    }
}


void Foam::diameterModels::populationBalanceModel::correctDilatationError()
{
    forAllIter
    (
        HashTable<volScalarField>,
        dilatationErrors_,
        iter
    )
    {
        volScalarField& dilatationError = iter();
        const word& phaseName = iter.key();
        const phaseModel& phase = fluid_.phases()[phaseName];
        const velocityGroup& velGroup = *velocityGroupPtrs_[phaseName];
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.thermo().rho();

        dilatationError =
            fvc::ddt(alpha) + fvc::div(phase.alphaPhi())
          - (fluid_.fvModels().source(alpha, rho) & rho)/rho;

        forAll(velGroup.sizeGroups(), i)
        {
            const sizeGroup& fi = velGroup.sizeGroups()[i];

            dilatationError -= Su_[fi.i()] - fvc::Sp(Sp_[fi.i()], fi);
        }
    }
}


void Foam::diameterModels::populationBalanceModel::calcAlphas()
{
    alphas_() = Zero;

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const phaseModel& phase = iter()->phase();

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

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const phaseModel& phase = iter()->phase();

        invDsm += max(phase, phase.residualAlpha())/(phase.d()*alphas_());
    }

    return 1/tInvDsm;
}


void Foam::diameterModels::populationBalanceModel::calcVelocity()
{
    U_() = Zero;

    forAllConstIter(HashTable<const velocityGroup*>, velocityGroupPtrs_, iter)
    {
        const phaseModel& phase = iter()->phase();

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
    dmdtfs_(),
    mesh_(fluid_.mesh()),
    name_(name),
    dict_
    (
        fluid_.subDict("populationBalanceCoeffs").subDict(name_)
    ),
    pimple_(mesh_.lookupObject<pimpleNoLoopControl>("solutionControl")),
    continuousPhase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict_.lookup("continuousPhase"))
        )
    ),
    sizeGroups_(),
    v_(),
    delta_(),
    Su_(),
    Sp_(),
    Sui_
    (
        IOobject
        (
            "Sui",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), Zero)
    ),
    coalescenceModels_
    (
        dict_.lookup("coalescenceModels"),
        coalescenceModel::iNew(*this)
    ),
    coalescenceRate_(),
    coalescencePairs_(),
    breakupModels_
    (
        dict_.lookup("breakupModels"),
        breakupModel::iNew(*this)
    ),
    breakupRate_(),
    binaryBreakupModels_
    (
        dict_.lookup("binaryBreakupModels"),
        binaryBreakupModel::iNew(*this)
    ),
    binaryBreakupRate_(),
    binaryBreakupPairs_(),
    drift_
    (
        dict_.lookup("driftModels"),
        driftModel::iNew(*this)
    ),
    driftRate_(),
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

    if (sizeGroups().size() < 3)
    {
        FatalErrorInFunction
            << "The populationBalance " << name_
            << " requires a minimum number of three sizeGroups to be specified."
            << exit(FatalError);
    }

    this->initialiseDmdtfs();

    if (coalescenceModels_.size() != 0)
    {
        coalescenceRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                     "coalescenceRate",
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
            new volScalarField
            (
                IOobject
                (
                    "breakupRate",
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
                    "binaryBreakupRate",
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

        calcDeltas();

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

    if (drift_.size() != 0)
    {
        driftRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "driftRate",
                    fluid_.time().name(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
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
                    fluid_.time().name(),
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

    if (velocityGroupPtrs_.size() > 1)
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
    const dimensionedScalar& x0 = sizeGroups()[0].x();
    const dimensionedScalar& xi = sizeGroups()[i].x();
    const dimensionedScalar& xm = sizeGroups().last().x();
    dimensionedScalar lowerBoundary(x0);
    dimensionedScalar upperBoundary(xm);

    if (i != 0) lowerBoundary = sizeGroups()[i-1].x();

    if (i != sizeGroups().size() - 1) upperBoundary = sizeGroups()[i+1].x();

    if ((i == 0 && v < x0) || (i == sizeGroups().size() - 1 && v > xm))
    {
        return v/xi;
    }
    else if (v < lowerBoundary || v > upperBoundary)
    {
        return 0;
    }
    else if (v.value() == xi.value())
    {
        return 1;
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


void Foam::diameterModels::populationBalanceModel::solve()
{
    if (!solveOnFinalIterOnly() || pimple_.finalIter())
    {
        const label nCorr = this->nCorr();
        const scalar tolerance =
            mesh_.solution().solverDict(name_).lookup<scalar>("tolerance");

        const bool updateSrc = updateSources();

        if (nCorr > 0 && updateSrc)
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

            if (updateSrc)
            {
                sources();
            }

            correctDilatationError();

            maxInitialResidual = 0;

            forAll(sizeGroups(), i)
            {
                sizeGroup& fi = sizeGroups_[i];
                const phaseModel& phase = fi.phase();
                const volScalarField& alpha = phase;
                const volScalarField& rho = phase.thermo().rho();
                const volScalarField& dilatationError =
                    dilatationErrors_[phase.name()];

                fvScalarMatrix sizeGroupEqn
                (
                    fvm::ddt(alpha, fi)
                  + fvm::div(phase.alphaPhi(), fi)
                  - fvm::Sp(dilatationError, fi)
                ==
                    fvc::Su(Su_[i], fi)
                  - fvm::Sp(Sp_[i], fi)
                  + fluid_.fvModels().source(alpha, rho, fi)/rho
                  - correction
                    (
                        fvm::Sp
                        (
                            max(phase.residualAlpha() - alpha, scalar(0))
                           /this->mesh().time().deltaT(),
                            fi
                        )
                    )
                );

                sizeGroupEqn.relax();
                fluid_.fvConstraints().constrain(sizeGroupEqn);

                maxInitialResidual = max
                (
                    sizeGroupEqn.solve().initialResidual(),
                    maxInitialResidual
                );

                fluid_.fvConstraints().constrain(fi);
            }
        }

        volScalarField fAlpha0
        (
            sizeGroups().first()*sizeGroups().first().phase()
        );

        volScalarField fAlphaN
        (
            sizeGroups().last()*sizeGroups().last().phase()
        );

        Info<< this->name() << " sizeGroup phase fraction first, last = "
            << fAlpha0.weightedAverage(this->mesh().V()).value()
            << ' ' << fAlphaN.weightedAverage(this->mesh().V()).value()
            << endl;
    }
}


void Foam::diameterModels::populationBalanceModel::correct()
{
    if (velocityGroupPtrs_.size() > 1)
    {
        calcAlphas();
        dsm_() = calcDsm();
        calcVelocity();
    }
}


// ************************************************************************* //
