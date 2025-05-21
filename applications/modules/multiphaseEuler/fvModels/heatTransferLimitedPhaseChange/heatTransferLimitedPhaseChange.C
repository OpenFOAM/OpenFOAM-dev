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

#include "heatTransferLimitedPhaseChange.H"
#include "fvmSup.H"
#include "multiphaseEuler.H"
#include "twoResistanceHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(heatTransferLimitedPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        heatTransferLimitedPhaseChange,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferLimitedPhaseChange::readCoeffs
(
    const dictionary& dict
)
{
    saturationModelPtr_.reset
    (
        saturationTemperatureModel::New
        (
            "saturationTemperature",
            dict
        ).ptr()
    );

    pressureImplicit_ = dict.lookup<bool>("pressureImplicit");

    if (pressureImplicit_ && !dmDotdpPtr_.valid())
    {
        dmDotdpPtr_.set
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    name() + ":dmDotdp",
                    mesh().time().name(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar(dimDensity/dimTime/dimPressure, 0)
            )
        );
    }

    if (!pressureImplicit_ && dmDotdpPtr_.valid())
    {
        dmDotdpPtr_.clear();
    }
}


void Foam::fv::heatTransferLimitedPhaseChange::correctMDot() const
{
    Info<< type() << ": " << name() << endl << incrIndent;

    const volScalarField::Internal& p = this->p();
    const volScalarField::Internal& T1 = phase1_.thermo().T();
    const volScalarField::Internal& T2 = phase2_.thermo().T();

    // Saturation temperature
    const volScalarField::Internal Tsat(saturationModelPtr_->Tsat(p));

    Info<< indent << "min/mean/max Tsat"
        << " = " << gMin(Tsat.primitiveField())
        << '/' << gAverage(Tsat.primitiveField())
        << '/' << gMax(Tsat.primitiveField())
        << endl;

    // Latent heat
    const volScalarField::Internal L(this->L(Tsat));

    // Heat transfer coefficients
    const Pair<tmp<volScalarField>> Hs =
        solver_.heatTransfer.Hs(phase1_, phase2_, scalar(0));
    const volScalarField::Internal& H1 = Hs.first();
    const volScalarField::Internal& H2 = Hs.second();

    // Relaxation factor
    const scalar f = mesh().solution().fieldRelaxationFactor(mDot_.member());

    // Recalculate the phase change rate
    mDot_ = (1 - f)*mDot_ + f*(H1*(T1 - Tsat) + H2*(T2 - Tsat))/L;

    Info<< indent << "min/mean/max mDot"
        << " = " << gMin(mDot_.primitiveField())
        << '/' << gAverage(mDot_.primitiveField())
        << '/' << gMax(mDot_.primitiveField())
        << endl;

    // Recalculate the phase change rate derivative w.r.t. pressure
    if (pressureImplicit_)
    {
        // Saturation temperature derivative w.r.t. pressure
        const volScalarField::Internal TsatPrime
        (
            saturationModelPtr_->TsatPrime(p)
        );

        dmDotdpPtr_() = (1 - f)*dmDotdpPtr_() - f*(H1 + H2)*TsatPrime/L;
    }

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferLimitedPhaseChange::heatTransferLimitedPhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange(name, modelType, mesh, dict, wordList()),
    solver_(mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName)),
    fluid_(solver_.fluid),
    phase1_(fluid_.phases()[phaseNames().first()]),
    phase2_(fluid_.phases()[phaseNames().second()]),
    saturationModelPtr_(nullptr),
    pressureImplicit_(false),
    pressureEquationIndex_(-1),
    mDot_
    (
        IOobject
        (
            name + ":mDot",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    dmDotdpPtr_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::heatTransferLimitedPhaseChange::Lfraction() const
{
    // Heat transfer coefficients
    const Pair<tmp<volScalarField>> Hs =
        solver_.heatTransfer.Hs(phase1_, phase2_);
    const volScalarField::Internal& H1 = Hs.first();
    const volScalarField::Internal& H2 = Hs.second();

    return H2/(H1 + H2);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::heatTransferLimitedPhaseChange::mDot() const
{
    return mDot_;
}


void Foam::fv::heatTransferLimitedPhaseChange::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &phase1_ || &alpha == &phase2_)
     && (&rho == &phase1_.rho() || &rho == &phase2_.rho())
     && &eqn.psi() == &solver_.p_rgh
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // the current phase loop
        if (pressureEquationIndex_ % 2 == 0) correctMDot();
        pressureEquationIndex_ ++;

        // Add linearisation if necessary
        if (pressureImplicit_)
        {
            const label s = &alpha == &phase1_ ? -1 : +1;

            eqn += correction(fvm::Sp(s*dmDotdpPtr_(), eqn.psi()));
        }
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


void Foam::fv::heatTransferLimitedPhaseChange::correct()
{
    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctMDot();
}


bool Foam::fv::heatTransferLimitedPhaseChange::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
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
