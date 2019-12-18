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

#include "greyMeanCombustion.H"
#include "combustionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(greyMeanCombustion, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        greyMeanCombustion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMeanCombustion::
greyMeanCombustion
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    greyMean(dict, mesh, typeName),
    EhrrCoeff_(coeffsDict_.lookup<scalar>("EhrrCoeff"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMeanCombustion::
~greyMeanCombustion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMeanCombustion::ECont
(
    const label bandI
) const
{
    tmp<volScalarField> E = greyMean::ECont(bandI);

    const word& name = combustionModel::combustionPropertiesName;
    E.ref() += EhrrCoeff_*mesh_.lookupObject<combustionModel>(name).Qdot();

    return E;
}


// ************************************************************************* //
