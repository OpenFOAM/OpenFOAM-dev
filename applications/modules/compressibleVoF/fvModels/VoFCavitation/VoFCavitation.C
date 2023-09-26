/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "VoFCavitation.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        namespace compressible
        {
            defineTypeNameAndDebug(VoFCavitation, 0);

            addToRunTimeSelectionTable
            (
                fvModel,
                VoFCavitation,
                dictionary
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::compressible::VoFCavitation::VoFCavitation
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),

    mixture_
    (
        mesh.lookupObjectRef<compressibleTwoPhaseVoFMixture>
        (
            "phaseProperties"
        )
    ),

    cavitation_(Foam::compressible::cavitationModel::New(dict, mixture_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::compressible::VoFCavitation::addSupFields() const
{
    return
    {
        mixture_.rho1().name(),
        mixture_.rho2().name()
    };
}


void Foam::fv::compressible::VoFCavitation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&rho == &mixture_.rho1() || &rho == &mixture_.rho2())
    {
        const volScalarField& alpha1 = mixture_.alpha1();
        const volScalarField& alpha2 = mixture_.alpha2();

        const scalar s = &alpha == &alpha1 ? +1 : -1;

        // Volume-fraction linearisation
        if (&alpha == &eqn.psi())
        {
            const Pair<tmp<volScalarField::Internal>> mDot12Alpha
            (
                cavitation_->mDot12Alpha()
            );
            const volScalarField::Internal& mDot1Alpha2 = mDot12Alpha[0]();
            const volScalarField::Internal& mDot2Alpha1 = mDot12Alpha[1]();

            eqn +=
                (&alpha == &alpha1 ? mDot1Alpha2 : mDot2Alpha1)
              - fvm::Sp(mDot1Alpha2 + mDot2Alpha1, eqn.psi());
        }

        // Pressure linearisation
        else if (eqn.psi().name() == "p_rgh")
        {
            const Pair<tmp<volScalarField::Internal>> mDot12P
            (
                cavitation_->mDot12P()
            );
            const volScalarField::Internal& mDot1P = mDot12P[0];
            const volScalarField::Internal& mDot2P = mDot12P[1];

            const volScalarField::Internal& rho =
                mesh().lookupObject<volScalarField>("rho");
            const volScalarField::Internal& gh =
                mesh().lookupObject<volScalarField>("gh");

            eqn +=
                fvm::Sp(s*(mDot1P - mDot2P), eqn.psi())
              + s*(mDot1P - mDot2P)*rho*gh
              - s*(mDot1P*cavitation_->pSat1() - mDot2P*cavitation_->pSat2());
        }

        // Explicit non-linearised value. Used in density predictors and
        // continuity error terms.
        else
        {
            const Pair<tmp<volScalarField::Internal>> mDot12Alpha
            (
                cavitation_->mDot12Alpha()
            );
            const volScalarField::Internal mDot1(mDot12Alpha[0]*alpha2);
            const volScalarField::Internal mDot2(mDot12Alpha[1]*alpha1);

            eqn += s*(mDot1 - mDot2);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << alpha.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::compressible::VoFCavitation::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::compressible::VoFCavitation::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::compressible::VoFCavitation::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::compressible::VoFCavitation::movePoints()
{
    return true;
}


void Foam::fv::compressible::VoFCavitation::correct()
{
    cavitation_->correct();
}


// ************************************************************************* //
