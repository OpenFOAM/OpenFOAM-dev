/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "twoResistanceHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"
#include "generateBlendedInterfacialModels.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(twoResistanceHeatTransfer, 0);
    addToRunTimeSelectionTable(fvModel, twoResistanceHeatTransfer, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::twoResistanceHeatTransfer::readCoeffs(const dictionary& dict)
{
    // Read and construct the models
    models_.clear();
    modelsTable tModels
    (
        generateBlendedInterfacialModels<blendedSidedHeatTransferModel>
        (
            fluid_,
            dict,
            fvModel::keywords,
            false
        )
    );
    models_.transfer(tModels);

    // Set the list of affected energy field names
    wordHashSet heNameSet;
    forAllConstIter(modelsTable, models_, modelIter)
    {
        forAllConstIter(phaseInterface, modelIter()->interface(), phaseIter)
        {
            heNameSet.insert(phaseIter().thermo().he().name());
        }
    }
    heNames_ = heNameSet.toc();
}


const Foam::fv::twoResistanceHeatTransfer::KsTable&
Foam::fv::twoResistanceHeatTransfer::Ks() const
{
    if (Ks_.empty())
    {
        forAllConstIter(modelsTable, models_, modelIter)
        {
            const phaseInterface& interface = modelIter()->interface();

            Pair<tmp<volScalarField>> interfaceKs
            (
                modelIter()->KinThe(interface.phase1()),
                modelIter()->KinThe(interface.phase2())
            );

            Ks_.insert(interface, interfaceKs);
        }
    }

    return Ks_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::twoResistanceHeatTransfer::twoResistanceHeatTransfer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    fluid_(mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)),
    models_(),
    heNames_(),
    energyEquationIndex_(-1),
    Ks_()
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::fv::twoResistanceHeatTransfer::Ks
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const phaseInterface interface(phase1, phase2);

    auto iter = Ks().find(interface);

    if (iter == Ks().end()) return Pair<tmp<volScalarField>>();

    Pair<tmp<volScalarField>> Ks(iter().first()(), iter().second()());

    return interface.index(phase2) ? Ks : reverse(Ks);
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::fv::twoResistanceHeatTransfer::Ks
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const scalar residualAlpha
) const
{
    const phaseInterface interface(phase1, phase2);

    auto iter = models_.find(interface);

    if (iter == models_.end()) return Pair<tmp<volScalarField>>();

    return
        Pair<tmp<volScalarField>>
        (
            iter()->KinThe(phase1, residualAlpha),
            iter()->KinThe(phase2, residualAlpha)
        );
}


Foam::wordList Foam::fv::twoResistanceHeatTransfer::addSupFields() const
{
    return heNames_;
}


void Foam::fv::twoResistanceHeatTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    // Ensure that the coefficients are up-to date if this is the first call in
    // the current phase loop
    if (energyEquationIndex_ % heNames_.size() == 0) Ks_.clear();
    energyEquationIndex_ ++;

    // Add sources to the energy equation of this phase
    const phaseModel& phase = fluid_.phases()[eqn.psi().group()];
    forAllConstIter(KsTable, Ks(), KIter)
    {
        const phaseInterface interface(fluid_, KIter.key());

        if (!interface.contains(phase)) continue;

        const phaseModel& otherPhase = interface.otherPhase(phase);

        const volScalarField& T = phase.thermo().T();
        const volScalarField& otherT = otherPhase.thermo().T();

        const volScalarField& K = KIter()[interface.index(phase)];
        const volScalarField& otherK = KIter()[interface.index(otherPhase)];

        const volScalarField KEff(K*otherK/(K + otherK));

        const volScalarField& Cpv = phase.thermo().Cpv();

        eqn += KEff*(otherT - T) + K/Cpv*he - fvm::Sp(K/Cpv, he);
    }
}


bool Foam::fv::twoResistanceHeatTransfer::movePoints()
{
    return true;
}


void Foam::fv::twoResistanceHeatTransfer::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::twoResistanceHeatTransfer::mapMesh(const polyMeshMap&)
{}


void Foam::fv::twoResistanceHeatTransfer::distribute(const polyDistributionMap&)
{}


void Foam::fv::twoResistanceHeatTransfer::correct()
{
    energyEquationIndex_ = 0;

    Ks_.clear();
}


bool Foam::fv::twoResistanceHeatTransfer::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
