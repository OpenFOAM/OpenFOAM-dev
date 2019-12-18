/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "wideBandCombustion.H"
#include "combustionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(wideBandCombustion, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        wideBandCombustion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wideBandCombustion::
wideBandCombustion
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    wideBand(dict, mesh, typeName)
{
    label bandi = 0;
    forAllConstIter(dictionary, coeffsDict_, iter)
    {
        if (!iter().isDict()) continue;

        iter().dict().lookup("EhrrCoeff") >> iEhrrCoeffs_[bandi];

        ++ bandi;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wideBandCombustion::
~wideBandCombustion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::wideBandCombustion::ECont
(
    const label bandi
) const
{
    tmp<volScalarField> E = wideBand::ECont(bandi);

    const word& name = combustionModel::combustionPropertiesName;
    E.ref() +=
        iEhrrCoeffs_[bandi]
       *mesh_.lookupObject<combustionModel>(name).Qdot()
       *(iBands_[bandi][1] - iBands_[bandi][0])
       /totalWaveLength_;

    return E;
}


// ************************************************************************* //
