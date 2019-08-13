/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
        popBal_.mesh().lookupObjectRef<reactingPhaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("reactingPhase"))
        )
    ),
    specie_
    (
        popBal_.mesh().lookupObjectRef<volScalarField>
        (
            IOobject::groupName
            (
                word(dict.lookup("specie")),
                reactingPhase_.name()
            )
        )
    ),
    nDmdt_
    (
        IOobject
        (
            "massChange",
            popBal.time().timeName(),
            popBal.mesh()
        ),
        popBal.mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    )
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

//- Add the corresponding species sources
void
Foam::diameterModels::nucleationModels::reactionDriven::registerPair
(
    populationBalanceModel::speciesDmdtTable& speciesDmdt
) const
{
    const phasePairKey key
    (
        velGroup_.phase().name(),
        reactingPhase_.name(),
        false
    );

    if (!speciesDmdt.found(key))
    {
        speciesDmdt.insert(key, new HashPtrTable<volScalarField>());
    }

    word specieName = IOobject::member(specie_.name());

    if (!speciesDmdt[key]->found(specieName))
    {
        speciesDmdt[key]->insert
        (
            specieName,
            new volScalarField
            (
                IOobject
                (
                    specieName+"NucleationSource",
                    popBal_.mesh().time().timeName(),
                    popBal_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                nDmdt_
            )
        );
    }
}


void Foam::diameterModels::nucleationModels::reactionDriven::correct()
{
    nDmdt_ = reactingPhase_*reactingPhase_.R(specie_) & specie_;
}


void
Foam::diameterModels::nucleationModels::reactionDriven::addToNucleationRate
(
    volScalarField& nucleationRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const volScalarField& rho = fi.phase().rho();

    nucleationRate +=
        popBal_.eta(i, pi/6.0*pow3(dNuc_))*nDmdt_/rho/fi.x();
}


void
Foam::diameterModels::nucleationModels::reactionDriven::addToSpeciesDmDt
(
    populationBalanceModel::speciesDmdtTable& speciesDmdt
) const
{
    const phasePairKey key
    (
        velGroup_.phase().name(),
        reactingPhase_.name(),
        false
    );

    word specieName = IOobject::member(specie_.name());

    *(*speciesDmdt[key])[specieName] += nDmdt_;
}


// ************************************************************************* //
