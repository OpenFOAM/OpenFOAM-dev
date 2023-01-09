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
#include "incompressibleTwoPhaseVoFMixture.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFCavitation::VoFCavitation
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(sourceName, modelType, dict, mesh),

    mixture_
    (
        mesh.lookupObjectRef<incompressibleTwoPhaseVoFMixture>
        (
            "phaseProperties"
        )
    ),

    cavitation_(cavitationModel::New(dict, mixture_)),

    alphaName_(mixture_.alpha1().name())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFCavitation::addSupFields() const
{
    return {alphaName_, "p_rgh", "U"};
}


void Foam::fv::VoFCavitation::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == alphaName_)
    {
        const volScalarField::Internal alpha1Coeff
        (
            1.0/mixture_.rho1()
          - mixture_.alpha1()()*(1.0/mixture_.rho1() - 1.0/mixture_.rho2())
        );

        const Pair<tmp<volScalarField::Internal>> mDot12Alpha
        (
            cavitation_->mDot12Alpha()
        );

        const volScalarField::Internal vDot1Alpha(alpha1Coeff*mDot12Alpha[0]());
        const volScalarField::Internal vDot2Alpha(alpha1Coeff*mDot12Alpha[1]());

        eqn += fvm::Sp(-vDot2Alpha - vDot1Alpha, eqn.psi()) + vDot1Alpha;
    }
}


void Foam::fv::VoFCavitation::addSup
(
    const volScalarField&,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == "p_rgh")
    {
        const volScalarField::Internal& rho =
            mesh().lookupObject<volScalarField>("rho");

        const volScalarField::Internal& gh =
            mesh().lookupObject<volScalarField>("gh");

        const dimensionedScalar pCoeff
        (
            1.0/mixture_.rho1() - 1.0/mixture_.rho2()
        );

        const Pair<tmp<volScalarField::Internal>> mDot12P
        (
            cavitation_->mDot12P()
        );

        const volScalarField::Internal vDot1P(pCoeff*mDot12P[0]);
        const volScalarField::Internal vDot2P(pCoeff*mDot12P[1]);

        eqn +=
            (vDot2P - vDot1P)*(cavitation_->pSat() - rho*gh)
          - fvm::Sp(vDot2P - vDot1P, eqn.psi());
    }
}


void Foam::fv::VoFCavitation::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == "U")
    {
        const surfaceScalarField& rhoPhi =
            mesh().lookupObject<surfaceScalarField>("rhoPhi");

        eqn += fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), eqn.psi());
    }
}


void Foam::fv::VoFCavitation::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::VoFCavitation::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::VoFCavitation::distribute(const polyDistributionMap&)
{}


bool Foam::fv::VoFCavitation::movePoints()
{
    return true;
}


void Foam::fv::VoFCavitation::correct()
{
    cavitation_->correct();
}


// ************************************************************************* //
