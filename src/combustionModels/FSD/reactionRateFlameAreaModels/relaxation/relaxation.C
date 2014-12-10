/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "relaxation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateFlameAreaModels
{
    defineTypeNameAndDebug(relaxation, 0);
    addToRunTimeSelectionTable
    (
        reactionRateFlameArea,
        relaxation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateFlameAreaModels::relaxation::relaxation
(
    const word modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
:
    reactionRateFlameArea(modelType, dict, mesh, combModel),
    correlation_(dict.subDict(typeName + "Coeffs").subDict(fuel_)),
    C_(readScalar(dict.subDict(typeName + "Coeffs").lookup("C"))),
    alpha_(readScalar(dict.subDict(typeName + "Coeffs").lookup("alpha")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateFlameAreaModels::relaxation::~relaxation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionRateFlameAreaModels::relaxation::correct
(
    const volScalarField& sigma
)
{

    dimensionedScalar omega0
    (
        "omega0",
        dimensionSet(1, -2, -1, 0, 0, 0, 0),
        correlation_.omega0()
    );

    dimensionedScalar sigmaExt
    (
        "sigmaExt",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        correlation_.sigmaExt()
    );

    dimensionedScalar omegaMin
    (
        "omegaMin",
        omega0.dimensions(),
        1e-4
    );

    const compressible::LESModel& lesModel =
        omega_.db().lookupObject<compressible::LESModel>("LESProperties");

    // Total strain : resolved and sub-grid (just LES for now)
    const volScalarField sigmaTotal
    (
        sigma + alpha_*lesModel.epsilon()/(lesModel.k() + lesModel.kMin())
    );

    const volScalarField omegaInf(correlation_.omega0Sigma(sigmaTotal));

    dimensionedScalar sigma0("sigma0", sigma.dimensions(), 0.0);

    const volScalarField tau(C_*mag(sigmaTotal));

    volScalarField Rc
    (
        (tau*omegaInf*(omega0 - omegaInf) + sqr(omegaMin)*sigmaExt)
       /(sqr(omega0 - omegaInf) + sqr(omegaMin))
    );

    const volScalarField rho(combModel_.rho());
    const surfaceScalarField phi(combModel_.phi());

    solve
    (
         fvm::ddt(rho, omega_)
       + fvm::div(phi, omega_, "div(phi,omega)")
      ==
         rho*Rc*omega0
       - fvm::SuSp(rho*(tau + Rc), omega_)
    );

    omega_.min(omega0);
    omega_.max(0.0);
}


bool  Foam::reactionRateFlameAreaModels::relaxation::read
(
    const dictionary& dict
)
{
    if (reactionRateFlameArea::read(dict))
    {
        coeffDict_ = dict.subDict(typeName + "Coeffs");
        coeffDict_.lookup("C") >> C_;
        coeffDict_.lookup("alpha") >> alpha_;
        correlation_.read
        (
            coeffDict_.subDict(fuel_)
        );
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
