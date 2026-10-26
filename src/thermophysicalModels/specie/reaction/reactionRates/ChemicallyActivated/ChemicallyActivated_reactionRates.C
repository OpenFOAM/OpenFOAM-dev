/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "makeRateReaction.H"

#include "ChemicallyActivated_reactionRate.H"
#include "Arrhenius_reactionRate.H"

#include "Lindemann_fallOffFunction.H"
#include "Troe_fallOffFunction.H"
#include "SRI_fallOffFunction.H"

#include "forGases.H"
#include "forLiquids.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace reactionRates
    {
        typedef ChemicallyActivated
        <
            Arrhenius,
            fallOffFunctions::Lindemann
        > ChemicallyActivatedArrheniusLindemann;
    }

    makeIrreversibleReaction(dictionary, ChemicallyActivatedArrheniusLindemann);
    forCoeffGases
    (
        makeReversibleReaction,
        dictionary,
        ChemicallyActivatedArrheniusLindemann
    );
    forCoeffLiquids
    (
        makeReversibleReaction,
        dictionary,
        ChemicallyActivatedArrheniusLindemann
    );

    namespace reactionRates
    {
        typedef ChemicallyActivated
        <
            Arrhenius,
            fallOffFunctions::Troe
        > ChemicallyActivatedArrheniusTroe;
    }

    makeIrreversibleReaction(dictionary, ChemicallyActivatedArrheniusTroe);
    forCoeffGases
    (
        makeReversibleReaction,
        dictionary,
        ChemicallyActivatedArrheniusTroe
    );
    forCoeffLiquids
    (
        makeReversibleReaction,
        dictionary,
        ChemicallyActivatedArrheniusTroe
    );

    namespace reactionRates
    {
        typedef ChemicallyActivated
        <
            Arrhenius,
            fallOffFunctions::SRI
        > ChemicallyActivatedArrheniusSRI;
    }

    makeIrreversibleReaction(dictionary, ChemicallyActivatedArrheniusSRI);
    forCoeffGases
    (
        makeReversibleReaction,
        dictionary,
        ChemicallyActivatedArrheniusSRI
    );
    forCoeffLiquids
    (
        makeReversibleReaction,
        dictionary,
        ChemicallyActivatedArrheniusSRI
    );
}

// ************************************************************************* //
