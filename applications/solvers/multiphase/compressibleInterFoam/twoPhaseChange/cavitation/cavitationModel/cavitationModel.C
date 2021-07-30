/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "cavitationModel.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeModels
{
    defineTypeNameAndDebug(cavitationModel, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::cavitationModel::cavitationModel
(
    const word& type,
    const compressibleTwoPhaseMixture& mixture
)
:
    twoPhaseChangeModel(type, mixture),
    pSat_("pSat", dimPressure, lookup("pSat"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::cavitationModel::Salpha
(
    volScalarField& alpha
) const
{
    const volScalarField::Internal alphalCoeff
    (
        1.0/rho1()
      - mixture_.alpha1()()*(1.0/rho1() - 1.0/rho2())
    );
    const Pair<tmp<volScalarField::Internal>> mDotAlphal = this->mDotAlphal();

    const volScalarField::Internal vDotcAlphal(alphalCoeff*mDotAlphal[0]());
    const volScalarField::Internal vDotvAlphal(alphalCoeff*mDotAlphal[1]());

    return Pair<tmp<volScalarField::Internal>>
    (
        1.0*vDotcAlphal,
        vDotvAlphal - vDotcAlphal
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::cavitationModel::Sp_rgh
(
    const volScalarField& rho,
    const volScalarField& gh,
    volScalarField& p_rgh
) const
{
    const volScalarField::Internal pCoeff(1.0/rho1() - 1.0/rho2());
    const Pair<tmp<volScalarField::Internal>> mDotP = this->mDotP();

    const volScalarField::Internal vDotcP(pCoeff*mDotP[0]);
    const volScalarField::Internal vDotvP(pCoeff*mDotP[1]);

    return
        (vDotvP - vDotcP)*(pSat() - rho()*gh())
      - fvm::Sp(vDotvP - vDotcP, p_rgh);
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::cavitationModel::SU
(
    const volScalarField& rho,
    const surfaceScalarField& rhoPhi,
    volVectorField& U
) const
{
    return fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U);
}


bool Foam::twoPhaseChangeModels::cavitationModel::read()
{
    if (twoPhaseChangeModel::read())
    {
        lookup("pSat") >> pSat_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
