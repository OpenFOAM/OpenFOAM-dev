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

#include "FallOff_reactionRate.H"
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
        typedef FallOff
        <
            Arrhenius,
            fallOffFunctions::Lindemann
        > FallOffArrheniusLindemann;
    }

    makeIrreversibleReaction(dictionary, FallOffArrheniusLindemann);
    forCoeffGases
    (
        makeReversibleReaction,
        dictionary,
        FallOffArrheniusLindemann
    );
    forCoeffLiquids
    (
        makeReversibleReaction,
        dictionary,
        FallOffArrheniusLindemann
    );

    namespace reactionRates
    {
        typedef FallOff<Arrhenius, fallOffFunctions::Troe> FallOffArrheniusTroe;
    }

    makeIrreversibleReaction(dictionary, FallOffArrheniusTroe);
    forCoeffGases(makeReversibleReaction, dictionary, FallOffArrheniusTroe);
    forCoeffLiquids(makeReversibleReaction, dictionary, FallOffArrheniusTroe);

    namespace reactionRates
    {
        typedef FallOff<Arrhenius, fallOffFunctions::SRI> FallOffArrheniusSRI;
    }

    makeIrreversibleReaction(dictionary, FallOffArrheniusSRI);
    forCoeffGases(makeReversibleReaction, dictionary, FallOffArrheniusSRI);
    forCoeffLiquids(makeReversibleReaction, dictionary, FallOffArrheniusSRI);
}

// ************************************************************************* //
