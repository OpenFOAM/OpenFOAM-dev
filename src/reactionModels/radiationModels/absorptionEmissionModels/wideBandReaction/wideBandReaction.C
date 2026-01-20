/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "wideBandReaction.H"
#include "reactionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(wideBandReaction, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        wideBandReaction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wideBandReaction::
wideBandReaction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    wideBand(dict, mesh, typeName)
{
    label bandi = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict()) continue;

        iter().dict().lookup("EhrrCoeff") >> iEhrrCoeffs_[bandi];

        ++ bandi;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wideBandReaction::
~wideBandReaction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::wideBandReaction::ECont
(
    const label bandi
) const
{
    tmp<volScalarField> E = wideBand::ECont(bandi);

    const word& name = reactionModel::reactionPropertiesName;
    E.ref() +=
        iEhrrCoeffs_[bandi]
       *mesh_.lookupObject<reactionModel>(name).Qdot()
       *(iBands_[bandi][1] - iBands_[bandi][0])
       /totalWaveLength_;

    return E;
}


// ************************************************************************* //
