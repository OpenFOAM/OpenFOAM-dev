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

#include "cloudAbsorptionEmission.H"
#include "parcelCloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(cloud, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        cloud,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::cloud::cloud
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(mesh),
    cloudNames_(dict.lookup("cloudNames"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::cloud::~cloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::cloud::aDisp(const label) const
{
    tmp<volScalarField> ta
    (
        volScalarField::New
        (
            "a",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );

    forAll(cloudNames_, i)
    {
        const parcelCloud& c =
            mesh_.objectRegistry::lookupObject<parcelCloud>(cloudNames_[i]);

        ta.ref() += c.ap();
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::cloud::eDisp
(
    const label bandI
) const
{
    tmp<volScalarField> te
    (
        volScalarField::New
        (
            "e",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::cloud::EDisp
(
    const label bandI
) const
{
    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "E",
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
        )
    );

    forAll(cloudNames_, i)
    {
        const parcelCloud& c =
            mesh_.objectRegistry::lookupObject<parcelCloud>(cloudNames_[i]);

        tE.ref() += c.Ep();
    }

    // Total emission is 4 times the projected emission
    return 4*tE;
}


// ************************************************************************* //
