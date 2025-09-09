/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "heatTransferSystem.H"

#include "fvmSup.H"

#include "heatTransferModel.H"
#include "generateBlendedInterfacialModels.H"

#include "twoResistanceHeatTransfer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatTransferSystem, 0);
}


const Foam::word Foam::heatTransferSystem::propertiesName
(
    "heatTransfer"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::heatTransferSystem::io(const phaseSystem& fluid) const
{
    typeIOobject<IOdictionary> result
    (
        propertiesName,
        fluid.mesh().time().constant(),
        fluid.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    // Check for and warn about heat transfer models being in the old
    // location in constant/phaseProperties
    if (!result.headerOk())
    {
        result.readOpt() = IOobject::NO_READ;

        if
        (
            fluid.found(modelName<heatTransferModel>())
         || !fluid.thermalPhases().empty()
        )
        {
            WarningInFunction
                << "Specifying a heat transfer model entry - "
                << modelName<heatTransferModel>() << " - in "
                << fluid.relativeObjectPath() << " is deprecated. The contents "
                << "of this entry should now be specified in "
                << result.relativeObjectPath() << "." << endl;
        }
    }
    else
    {
        if (fluid.found(modelName<heatTransferModel>()))
        {
            WarningInFunction
                << "Heat transfer model entry - "
                << modelName<heatTransferModel>() << " - in "
                << fluid.relativeObjectPath() << " is no longer used. The "
                << "contents of this entry are now read from "
                << result.relativeObjectPath() << "." << endl;
        }
    }

    return result;
}


const Foam::dictionary& Foam::heatTransferSystem::modelsDict() const
{
    const word key = modelName<heatTransferModel>();

    return
        readOpt() == IOobject::NO_READ && fluid_.found(key)
      ? fluid_.subDict(key)
      : *this;
}


template<class ... Args>
Foam::Pair<Foam::tmp<Foam::volScalarField>> Foam::heatTransferSystem::Hs
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    Args ... args
) const
{
    auto error = [&](const word& why)
    {
        FatalErrorInFunction
            << why << " two-resistance heat transfer models found that "
            << "provide a heat transfer coefficient between phases "
            << phase1.name() << " and " << phase2.name() << exit(FatalError);
    };

    autoPtr<Pair<tmp<volScalarField>>> HsPtr;

    const phaseInterface interface(phase1, phase2);

    auto iter = sidedModels_.find(interface);

    if (iter != sidedModels_.end())
    {
        HsPtr.set
        (
            new Pair<tmp<volScalarField>>
            (
                iter()->KinThe(phase1, args ...),
                iter()->KinThe(phase2, args ...)
            )
        );
    }

    const Foam::fvModels& fvModels = Foam::fvModels::New(fluid_.mesh());

    forAll(fvModels, i)
    {
        if (!isA<fv::twoResistanceHeatTransfer>(fvModels[i])) continue;

        const fv::twoResistanceHeatTransfer& heatTransferFvModel =
            refCast<const fv::twoResistanceHeatTransfer>(fvModels[i]);

        Pair<tmp<volScalarField>> Hs =
            heatTransferFvModel.Ks(phase1, phase2, args ...);

        if (!Hs.first().valid() || !Hs.second().valid()) continue;

        if (HsPtr.valid()) error("Multiple");

        HsPtr.set(new Pair<tmp<volScalarField>>(Hs.first(), Hs.second()));
    }

    if (!HsPtr.valid()) error("No");

    return HsPtr();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferSystem::heatTransferSystem
(
    const phaseSystem& fluid
)
:
    IOdictionary(io(fluid)),
    fluid_(fluid),
    models_(),
    sidedModels_()
{
    // If we have no entries and there are no thermal phases then this case
    // does not need any heat transfer modelling and we can just quit without
    // trying to construct anything
    if
    (
        readOpt() == IOobject::NO_READ
     && !fluid_.found(modelName<heatTransferModel>())
     && fluid.thermalPhases().empty()
    ) return;

    modelsTable models
    (
        generateBlendedInterfacialModels<blendedHeatTransferModel>
        (
            fluid,
            modelsDict(),
            wordHashSet(),
            true
        )
    );

    models_.transfer(models);

    sidedModelsTable sidedModels
    (
        generateBlendedInterfacialModels<blendedSidedHeatTransferModel>
        (
            fluid,
            modelsDict(),
            wordHashSet(),
            true
        )
    );

    sidedModels_.transfer(sidedModels);

    forAllConstIter(modelsTable, models_, modelIter)
    {
        if (sidedModels_.found(modelIter.key()))
        {
            const phaseInterface interface(fluid_, modelIter.key());

            FatalIOErrorInFunction(modelsDict())
                << "One-resistance and two-resistance heat transfer models "
                << "both specified between phases "
                << interface.phase1().name() << " and "
                << interface.phase2().name() << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferSystem::~heatTransferSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>> Foam::heatTransferSystem::Hs
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return Hs<>(phase1, phase2);
}


Foam::Pair<Foam::tmp<Foam::volScalarField>> Foam::heatTransferSystem::Hs
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const scalar residualAlpha
) const
{
    return Hs<scalar>(phase1, phase2, residualAlpha);
}


Foam::autoPtr<Foam::HashPtrTable<Foam::fvScalarMatrix>>
Foam::heatTransferSystem::heatTransfer() const
{
    autoPtr<HashPtrTable<fvScalarMatrix>> eqnsPtr
    (
        new HashPtrTable<fvScalarMatrix>()
    );
    HashPtrTable<fvScalarMatrix>& eqns = eqnsPtr();

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    forAllConstIter(modelsTable, models_, modelIter)
    {
        const phaseInterface interface(fluid_, modelIter.key());

        const volScalarField H(modelIter()->K());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField& he = phase.thermo().he();
            const volScalarField Cpv(phase.thermo().Cpv());

            const volScalarField Hstabilised
            (
                iter.otherPhase()
               /max(iter.otherPhase(), iter.otherPhase().residualAlpha())
               *H
            );

            *eqns[phase.name()] +=
                Hstabilised*(otherPhase.thermo().T() - phase.thermo().T())
              + Hstabilised/Cpv*he - fvm::Sp(Hstabilised/Cpv, he);
        }
    }

    forAllConstIter(sidedModelsTable, sidedModels_, sidedModelIter)
    {
        const phaseInterface interface(fluid_, sidedModelIter.key());

        Pair<volScalarField> Hs
        (
            sidedModelIter()->KinThe(interface.phase1()),
            sidedModelIter()->KinThe(interface.phase2())
        );

        const volScalarField HEff
        (
            Hs.first()*Hs.second()/(Hs.first() + Hs.second())
        );

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField& he = phase.thermo().he();
            const volScalarField Cpv(phase.thermo().Cpv());

            const volScalarField& H = Hs[iter.index()];

            *eqns[phase.name()] +=
                HEff*(otherPhase.thermo().T() - phase.thermo().T())
              + H/Cpv*he - fvm::Sp(H/Cpv, he);
        }
    }

    return eqnsPtr;
}


bool Foam::heatTransferSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
