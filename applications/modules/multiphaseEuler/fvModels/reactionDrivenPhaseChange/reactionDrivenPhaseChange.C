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

#include "reactionDrivenPhaseChange.H"
#include "multicomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(reactionDrivenPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        reactionDrivenPhaseChange,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::reactionDrivenPhaseChange::readCoeffs(const dictionary& dict)
{
    reReadSpecies(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::reactionDrivenPhaseChange::reactionDrivenPhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange
    (
        name,
        modelType,
        mesh,
        dict,
        readSpecies(coeffs(modelType, dict), true)
    ),
    fluid_(mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)),
    phase1_(fluid_.phases()[phaseNames().first()]),
    phase2_(fluid_.phases()[phaseNames().second()])
{
    readCoeffs(coeffs(dict));

    const ThermoRefPair<multicomponentThermo> mcThermos =
        thermos().thermos<multicomponentThermo>();

    if (!mcThermos.either())
    {
        WarningInFunction
            << "Model " << name << " of type " << modelType
            << " applied to two pure phases " << phase1_.name()
            << " and " << phase2_.name() << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::reactionDrivenPhaseChange::mDot(const label mDoti) const
{
    const volScalarField::Internal& alpha1 = phase1_, alpha2 = phase2_;

    const labelPair specieis = this->specieis(mDoti);

    const ThermoRefPair<multicomponentThermo> mcThermos =
        thermos().thermos<multicomponentThermo>();

    tmp<volScalarField::Internal> tResult =
        volScalarField::Internal::New
        (
            name() + ":mDot_" + species()[mDoti],
            mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );

    if (mcThermos.valid().first())
    {
        tResult.ref() += alpha1*phase1_.R(specieis.first());
    }

    if (mcThermos.valid().second())
    {
        tResult.ref() -= alpha2*phase2_.R(specieis.second());
    }

    return tResult;
}


bool Foam::fv::reactionDrivenPhaseChange::read(const dictionary& dict)
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
