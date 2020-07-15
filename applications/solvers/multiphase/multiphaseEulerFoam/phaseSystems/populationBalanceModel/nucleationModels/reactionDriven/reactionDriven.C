/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvmDdt.H"
#include "fvcDdt.H"

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
    velGroup_
    (
        refCast<const velocityGroup>
        (
            popBal.mesh().lookupObject<phaseModel>
            (
                IOobject::groupName
                (
                    "alpha",
                    dict.lookup("velocityGroup")
                )
            ).dPtr()()
        )
    ),
    reactingPhase_
    (
        popBal_.mesh().lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("reactingPhase"))
        )
    ),
    pair_
    (
        popBal_.fluid().phasePairs()
        [
            phasePair(velGroup_.phase(), reactingPhase_)
        ]
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

    const volScalarField& dmidtf =
        popBal_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                IOobject::groupName
                (
                    dmdtfName_,
                    specieName_
                ),
                pair_.name()
            )
        );

    const scalar dmidtfSign =
        velGroup_.phase().name() == pair_.first() ? +1 : -1;

    nucleationRate +=
        popBal_.eta(i, pi/6.0*pow3(dNuc_))*dmidtfSign*dmidtf/rho/fi.x();
}


// ************************************************************************* //
