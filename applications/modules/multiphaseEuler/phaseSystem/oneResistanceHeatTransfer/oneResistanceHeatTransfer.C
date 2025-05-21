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

#include "oneResistanceHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"
#include "generateBlendedInterfacialModels.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(oneResistanceHeatTransfer, 0);
    addToRunTimeSelectionTable(fvModel, oneResistanceHeatTransfer, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::oneResistanceHeatTransfer::readCoeffs(const dictionary& dict)
{
    // Read and construct the models
    models_.clear();
    modelsTable tModels
    (
        generateBlendedInterfacialModels<blendedHeatTransferModel>
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


const Foam::fv::oneResistanceHeatTransfer::KsTable&
Foam::fv::oneResistanceHeatTransfer::Ks() const
{
    if (Ks_.empty())
    {
        forAllConstIter(modelsTable, models_, modelIter)
        {
            Ks_.insert
            (
                modelIter()->interface(),
                modelIter()->K()
            );
        }
    }

    return Ks_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::oneResistanceHeatTransfer::oneResistanceHeatTransfer
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

Foam::tmp<Foam::volScalarField> Foam::fv::oneResistanceHeatTransfer::K
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    auto iter = Ks().find(phaseInterface(phase1, phase2));

    return
        iter == Ks().end()
      ? tmp<volScalarField>()
      : tmp<volScalarField>((*iter)());
}


Foam::wordList Foam::fv::oneResistanceHeatTransfer::addSupFields() const
{
    return heNames_;
}


void Foam::fv::oneResistanceHeatTransfer::addSup
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

        const volScalarField& K = KIter();

        const volScalarField stabilisedK
        (
            otherPhase/max(otherPhase, otherPhase.residualAlpha())*K
        );

        const volScalarField& Cpv = phase.thermo().Cpv();

        eqn += stabilisedK*(otherT - T + he/Cpv) - fvm::Sp(stabilisedK/Cpv, he);
    }
}


bool Foam::fv::oneResistanceHeatTransfer::movePoints()
{
    return true;
}


void Foam::fv::oneResistanceHeatTransfer::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::oneResistanceHeatTransfer::mapMesh(const polyMeshMap&)
{}


void Foam::fv::oneResistanceHeatTransfer::distribute(const polyDistributionMap&)
{}


void Foam::fv::oneResistanceHeatTransfer::correct()
{
    energyEquationIndex_ = 0;

    Ks_.clear();
}


bool Foam::fv::oneResistanceHeatTransfer::read(const dictionary& dict)
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
