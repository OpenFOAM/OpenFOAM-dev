/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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

#include "slurry.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(slurry, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        slurry,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::slurry::slurry
(
    const fvMesh& mesh,
    const word& group
)
:
    mixtureViscosityModel(mesh, group),
    alpha_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                lookupOrDefault<word>("alpha", "alpha"),
                group
            )
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::slurry::mu
(
    const volScalarField& muc,
    const volVectorField& U
) const
{
    return
    (
        muc*(1.0 + 2.5*alpha_ + 10.05*sqr(alpha_) + 0.00273*exp(16.6*alpha_))
    );
}


bool Foam::mixtureViscosityModels::slurry::read()
{
    return mixtureViscosityModel::read();
}


// ************************************************************************* //
