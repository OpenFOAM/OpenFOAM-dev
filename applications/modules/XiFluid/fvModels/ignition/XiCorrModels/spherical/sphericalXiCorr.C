/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "sphericalXiCorr.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiCorrModels
{
    defineTypeNameAndDebug(spherical, 0);
    addToRunTimeSelectionTable(XiCorrModel, spherical, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiCorrModels::spherical::readCoeffs(const dictionary& dict)
{
    XiCorrModel::readCoeffs(dict);
    sphereFraction_.readIfPresent(dict);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiCorrModels::spherical::spherical
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    XiCorrModel(mesh, dict),
    sphereFraction_("sphereFraction", dimless, 1)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiCorrModels::spherical::~spherical()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::XiCorrModels::spherical::Ak
(
    const dimensionedScalar& Vk
) const
{
    // Radius of the ignition kernel
    dimensionedScalar rk
    (
        pow
        (
            (3.0/4.0)*Vk
           /(sphereFraction_*constant::mathematical::pi),
            1.0/3.0
        )
    );

    // Return area of the ignition kernel
    return sphereFraction_*4*constant::mathematical::pi*sqr(rk);
}


// ************************************************************************* //
