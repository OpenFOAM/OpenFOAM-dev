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
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseTransferModels
{
    defineTypeNameAndDebug(reactionDriven, 0);
    addToRunTimeSelectionTable(phaseTransferModel, reactionDriven, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<bool Index>
Foam::word Foam::phaseTransferModels::reactionDriven::speciesKey() const
{
    return
        IOobject::groupName
        (
            "species",
            (!Index ? interface_.phase1() : interface_.phase2()).name()
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransferModels::reactionDriven::reactionDriven
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    phaseTransferModel(dict, interface),
    interface_(interface),
    species1_(dict.lookupOrDefault<wordList>(speciesKey<0>(), wordList())),
    species2_(dict.lookupOrDefault<wordList>(speciesKey<1>(), wordList())),
    species_()
{
    if (!dict.found(speciesKey<0>()) && !dict.found(speciesKey<1>()))
    {
        FatalIOErrorInFunction(dict)
            << "No transferring species specified. Specify either "
            << speciesKey<0>() << " or " << speciesKey<1>() << " or both."
            << exit(FatalIOError);
    }

    wordList species(species1_);
    species.append(species2_);
    species_.transfer(species);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransferModels::reactionDriven::~reactionDriven()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::hashedWordList&
Foam::phaseTransferModels::reactionDriven::species() const
{
    return species_;
}


Foam::HashPtrTable<Foam::volScalarField>
Foam::phaseTransferModels::reactionDriven::dmidtf() const
{
    HashPtrTable<volScalarField> result;

    const phaseModel& phase1 = interface_.phase1();
    const phaseModel& phase2 = interface_.phase2();

    forAll(species1_, i)
    {
        volScalarField& Y1 =
            const_cast<volScalarField&>(phase1.Y(species1_[i]));
        result.set(species1_[i], (- phase1*phase1.R(Y1) & Y1).ptr());
    }

    forAll(species2_, i)
    {
        volScalarField& Y2 =
            const_cast<volScalarField&>(phase2.Y(species2_[i]));
        result.set(species2_[i], (phase2*phase2.R(Y2) & Y2).ptr());
    }

    return result;
};


// ************************************************************************* //
