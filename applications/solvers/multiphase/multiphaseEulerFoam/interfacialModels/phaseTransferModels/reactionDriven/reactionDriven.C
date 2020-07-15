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
#include "phasePair.H"
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransferModels::reactionDriven::reactionDriven
(
    const dictionary& dict,
    const phasePair& pair
)
:
    phaseTransferModel(dict, pair),
    reactingName_(dict.lookup("reactingPhase")),
    reactingPhase_
    (
        reactingName_ == pair_.first() ? pair_.phase1() : pair_.phase2()
    ),
    otherPhase_
    (
        pair.otherPhase(reactingPhase_)
    ),
    sign_
    (
        reactingName_ == pair_.first() ? -1 : 1
    ),
    species_(dict.lookup("species"))
{}


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

    forAll(species_, i)
    {
        const word name = species_[i];

        volScalarField& Y =
            const_cast<volScalarField&>(reactingPhase_.Y(name));

        result.set
        (
            species_[i],
            (sign_*reactingPhase_*reactingPhase_.R(Y) & Y).ptr()
        );
    }

    return result;
};


// ************************************************************************* //
