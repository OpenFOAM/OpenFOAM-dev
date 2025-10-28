/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "cylindricalKernelShape.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kernelShapes
{
    defineTypeNameAndDebug(cylindrical, 0);
    addToRunTimeSelectionTable(kernelShape, cylindrical, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::kernelShapes::cylindrical::readCoeffs(const dictionary& dict)
{
    kernelShape::readCoeffs(dict);
    thickness_.read(dict);
    cylinderFraction_.readIfPresent(dict);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kernelShapes::cylindrical::cylindrical
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    kernelShape(mesh, dict),
    thickness_("thickness", dimLength, 0),
    cylinderFraction_("cylinderFraction", dimless, 1)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kernelShapes::cylindrical::~cylindrical()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::kernelShapes::cylindrical::Ak
(
    const dimensionedScalar& Vk
) const
{
    // Radius of the ignition kernel
    dimensionedScalar rk
    (
        sqrt(Vk/(cylinderFraction_*constant::mathematical::pi*thickness_))
    );

    // Return area of the ignition kernel
    return cylinderFraction_*2*constant::mathematical::pi*rk*thickness_;
}


// ************************************************************************* //
