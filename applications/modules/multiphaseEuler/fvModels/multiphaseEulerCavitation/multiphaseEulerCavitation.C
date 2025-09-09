/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "multiphaseEulerCavitation.H"
#include "fvmSup.H"
#include "multiphaseEuler.H"
#include "cavitationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseEulerCavitation, 0);
    addToRunTimeSelectionTable(fvModel, multiphaseEulerCavitation, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::multiphaseEulerCavitation::readCoeffs
(
    const dictionary& dict
)
{
    cavitationModel_.reset
    (
        compressible::cavitationModel::New
        (
            dict,
            interface_,
            interface_.index(liquid_)
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseEulerCavitation::multiphaseEulerCavitation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange(name, modelType, mesh, dict, wordList()),
    fluid_(mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)),
    liquid_(fluid_.phases()[phaseNames().first()]),
    vapour_(fluid_.phases()[phaseNames().second()]),
    interface_(liquid_, vapour_),
    p_rgh_
    (
        mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName).p_rgh
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::multiphaseEulerCavitation::Lfraction() const
{
    // Put all the latent heat into the liquid
    return
        volScalarField::Internal::New
        (
            name() + ":Lfraction",
            mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::multiphaseEulerCavitation::mDot() const
{
    const Pair<tmp<volScalarField::Internal>> coeffs
    (
        cavitationModel_->mDot12P()
    );

    const volScalarField::Internal& p = liquid_.fluidThermo().p();
    const volScalarField::Internal pSat1(cavitationModel_->pSat1());
    const volScalarField::Internal pSat2(cavitationModel_->pSat2());

    return
        volScalarField::Internal::New
        (
            name() + ":mDot",
            coeffs[0]*(p - pSat1) - coeffs[1]*(p - pSat2)
        );
}


void Foam::fv::multiphaseEulerCavitation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &liquid_ || &alpha == &vapour_)
     && (&rho == &liquid_.rho() || &rho == &vapour_.rho())
     && &eqn.psi() == &p_rgh_
    )
    {
        const Pair<tmp<volScalarField::Internal>> coeffs
        (
            cavitationModel_->mDot12P()
        );

        eqn += correction(fvm::Sp(coeffs[0] - coeffs[1], eqn.psi()));
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


bool Foam::fv::multiphaseEulerCavitation::read(const dictionary& dict)
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
