/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "reactionDriven.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(reactionDriven, 0);
    addToRunTimeSelectionTable
    (
        nucleationModel,
        reactionDriven,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::nucleationModels::reactionDriven::
reactionDriven
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    nucleationModel(popBal, dict),
    dNuc_("nucleationDiameter", dimLength, dict),
    reactingPhase_
    (
        popBal_.fluid().phases()[dict.lookup<word>("reactingPhase")]
    ),
    dmdtfName_(dict.lookup("dmdtf")),
    specieName_(dict.lookup("specie"))
{
    if
    (
        dNuc_.value() < velGroup_.sizeGroups().first().dSph().value()
     || dNuc_.value() > velGroup_.sizeGroups().last().dSph().value()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Nucleation diameter " << dNuc_.value() << "m outside of range ["
            << velGroup_.sizeGroups().first().dSph().value() << ", "
            << velGroup_.sizeGroups().last().dSph().value() << "]." << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::nucleationModels::reactionDriven::addToNucleationRate
(
    volScalarField& nucleationRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const volScalarField& rho = fi.phase().rho();

    const phaseInterface interface(velGroup_.phase(), reactingPhase_);

    const volScalarField& dmidtf =
        popBal_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                IOobject::groupName(dmdtfName_, specieName_),
                interface.name()
            )
        );

    const scalar dmidtfSign =
        interface.index(velGroup_.phase()) == 0 ? +1 : -1;

    nucleationRate +=
        popBal_.eta(i, pi/6*pow3(dNuc_))*dmidtfSign*dmidtf/rho/fi.x();
}


// ************************************************************************* //
