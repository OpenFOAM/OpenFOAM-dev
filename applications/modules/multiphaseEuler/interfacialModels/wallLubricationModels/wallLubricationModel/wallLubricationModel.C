/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "wallLubricationModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallLubricationModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(wallLubricationModel, 0);
    defineRunTimeSelectionTable(wallLubricationModel, dictionary);
}

const Foam::dimensionSet Foam::wallLubricationModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModel::zeroGradWalls
(
    tmp<volVectorField> tFi
) const
{
    volVectorField& Fi = tFi.ref();
    const fvPatchList& patches =  Fi.mesh().boundary();

    volVectorField::Boundary& FiBf = Fi.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            fvPatchVectorField& Fiw = FiBf[patchi];
            Fiw = Fiw.patchInternalField();
        }
    }

    return tFi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModel::wallLubricationModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    wallDependentModel(interface.mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModel::~wallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::blendedWallLubricationModel::F() const
{
    return
        evaluate
        (
            &wallLubricationModel::F,
            "F",
            wallLubricationModel::dimF,
            true
        );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::blendedWallLubricationModel::Ff() const
{
    return
        evaluate
        (
            &wallLubricationModel::Ff,
            "Ff",
            wallLubricationModel::dimF*dimArea,
            true
        );
}


// ************************************************************************* //
