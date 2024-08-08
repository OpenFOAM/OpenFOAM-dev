/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(absorptionEmissionModel, 0);
        defineRunTimeSelectionTable(absorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModel::absorptionEmissionModel
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModel::~absorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::a(const label bandI) const
{
    return aDisp(bandI) + aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::aCont(const label bandI) const
{
    return volScalarField::New
    (
        "aCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::aDisp(const label bandI) const
{
    return volScalarField::New
    (
        "aDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::e(const label bandI) const
{
    return eDisp(bandI) + eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::eCont(const label bandI) const
{
    return volScalarField::New
    (
        "eCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::eDisp(const label bandI) const
{
    return volScalarField::New
    (
        "eDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::E(const label bandI) const
{
    return EDisp(bandI) + ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::ECont(const label bandI) const
{
    return volScalarField::New
    (
        "ECont",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModel::EDisp(const label bandI) const
{
    return volScalarField::New
    (
        "EDisp",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


Foam::label Foam::radiationModels::absorptionEmissionModel::nBands() const
{
    return pTraits<label>::one;
}


const Foam::Vector2D<Foam::scalar>&
Foam::radiationModels::absorptionEmissionModel::bands(const label n) const
{
    return Vector2D<scalar>::one;
}


bool Foam::radiationModels::absorptionEmissionModel::isGrey() const
{
    return false;
}


void Foam::radiationModels::absorptionEmissionModel::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
    a = this->a();
    aj[0] =  a;
}


// ************************************************************************* //
