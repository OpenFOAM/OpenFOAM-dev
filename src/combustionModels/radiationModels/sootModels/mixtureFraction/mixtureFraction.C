/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "mixtureFraction.H"
#include "singleStepCombustion.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{
    defineTypeNameAndDebug(mixtureFraction, 0);

    addToRunTimeSelectionTable
    (
        sootModel,
        mixtureFraction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::sootModels::mixtureFraction::mixtureFraction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),
    soot_
    (
        IOobject
        (
            "soot",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),
    nuSoot_(coeffsDict_.lookup<scalar>("nuSoot")),
    Wsoot_(coeffsDict_.lookup<scalar>("Wsoot")),
    sootMax_(-1),
    mappingFieldName_
    (
        coeffsDict_.lookupOrDefault<word>("mappingField", "none")
    ),
    mapFieldMax_(1)
{
    const combustionModels::singleStepCombustion& combustion =
        mesh.lookupObject<combustionModels::singleStepCombustion>
        (
            combustionModel::combustionPropertiesName
        );

    const reaction& singleReaction = combustion.singleReaction();

    scalar totalMol = 0;
    forAll(singleReaction.rhs(), i)
    {
        const scalar stoichCoeff = singleReaction.rhs()[i].stoichCoeff;
        totalMol += mag(stoichCoeff);
    }

    totalMol += nuSoot_;

    scalarList Xi(singleReaction.rhs().size());

    scalar Wm = 0;
    forAll(singleReaction.rhs(), i)
    {
        const label speciei = singleReaction.rhs()[i].index;
        const scalar stoichCoeff = singleReaction.rhs()[i].stoichCoeff;
        Xi[i] = mag(stoichCoeff)/totalMol;
        Wm += Xi[i]*combustion.thermo().WiValue(speciei);
    }

    scalarList Yprod0(combustion.thermo().species().size(), 0.0);

    forAll(singleReaction.rhs(), i)
    {
        const label speciei = singleReaction.rhs()[i].index;
        Yprod0[speciei] = combustion.thermo().WiValue(speciei)/Wm*Xi[i];
    }

    const scalar XSoot = nuSoot_/totalMol;
    Wm += XSoot*Wsoot_;

    sootMax_ = XSoot*Wsoot_/Wm;

    Info << "Maximum soot mass concentrations: " << sootMax_ << nl;

    if (mappingFieldName_ == "none")
    {
        const label index = singleReaction.rhs()[0].index;
        mappingFieldName_ = combustion.thermo().Y(index).name();
    }

    const label mapFieldIndex =
        combustion.thermo().species()[mappingFieldName_];

    mapFieldMax_ = Yprod0[mapFieldIndex];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::sootModels::mixtureFraction::~mixtureFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::sootModels::mixtureFraction::correct()
{
    const volScalarField& mapField =
        mesh_.lookupObject<volScalarField>(mappingFieldName_);

    soot_ = sootMax_*(mapField/mapFieldMax_);
}


// ************************************************************************* //
