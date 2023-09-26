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
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),

    mixture_
    (
        mesh.lookupObjectRef<incompressibleTwoPhaseVoFMixture>
        (
            "phaseProperties"
        )
    ),

    cavitation_(cavitationModel::New(dict, mixture_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFCavitation::addSupFields() const
{
    return
    {
        mixture_.alpha1().name(),
        mixture_.alpha2().name(),
        "U"
    };
}


void Foam::fv::VoFCavitation::addSup
(
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&alpha == &mixture_.alpha1() || &alpha == &mixture_.alpha2())
    {
        const volScalarField& alpha1 = mixture_.alpha1();
        const volScalarField& alpha2 = mixture_.alpha2();

        const dimensionedScalar& rho =
            &alpha == &alpha1 ? mixture_.rho1() : mixture_.rho2();

        const scalar s = &alpha == &alpha1 ? +1 : -1;

        // Volume-fraction linearisation
        if (&alpha == &eqn.psi())
        {
            const Pair<tmp<volScalarField::Internal>> mDot12Alpha
            (
                cavitation_->mDot12Alpha()
            );
            const volScalarField::Internal vDot1Alpha2(mDot12Alpha[0]/rho);
            const volScalarField::Internal vDot2Alpha1(mDot12Alpha[1]/rho);

            eqn +=
                (&alpha == &alpha1 ? vDot1Alpha2 : vDot2Alpha1)
              - fvm::Sp(vDot1Alpha2 + vDot2Alpha1, eqn.psi());
        }

        // Pressure linearisation
        else if (eqn.psi().name() == "p_rgh")
        {
            const Pair<tmp<volScalarField::Internal>> mDot12P
            (
                cavitation_->mDot12P()
            );
            const volScalarField::Internal vDot1P(mDot12P[0]/rho);
            const volScalarField::Internal vDot2P(mDot12P[1]/rho);

            const volScalarField::Internal& rho =
                mesh().lookupObject<volScalarField>("rho");
            const volScalarField::Internal& gh =
                mesh().lookupObject<volScalarField>("gh");

            eqn +=
                fvm::Sp(s*(vDot1P - vDot2P), eqn.psi())
              + s*(vDot1P - vDot2P)*rho*gh
              - s*(vDot1P - vDot2P)*cavitation_->pSat();
        }

        // Explicit non-linearised value. Probably not used.
        else
        {
            const Pair<tmp<volScalarField::Internal>> mDot12Alpha
            (
                cavitation_->mDot12Alpha()
            );
            const volScalarField::Internal vDot1(mDot12Alpha[0]*alpha2/rho);
            const volScalarField::Internal vDot2(mDot12Alpha[1]*alpha1/rho);

            eqn += s*(vDot1 - vDot2);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << alpha.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFCavitation::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (U.name() == "U")
    {
        const surfaceScalarField& rhoPhi =
            mesh().lookupObject<surfaceScalarField>("rhoPhi");

        if (&U == &eqn.psi())
        {
            eqn += fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), eqn.psi());
        }
        else
        {
            eqn += (fvc::ddt(rho) + fvc::div(rhoPhi))*U;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << U.name() << " is not implemented"
            << exit(FatalError);
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
