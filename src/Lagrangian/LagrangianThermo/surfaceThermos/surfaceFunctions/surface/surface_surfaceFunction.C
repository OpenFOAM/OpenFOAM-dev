/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "surface_surfaceFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFunctions
{
    defineTypeNameAndDebug(surface, 0);
    addToRunTimeSelectionTable(surfaceFunction, surface, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceFunctions::surface::thermalConductivityCoeffs
(
    const LagrangianSubMesh& subMesh,
    autoPtr<dimensionedScalar>& fValuePtr,
    tmp<LagrangianSubScalarField>& tfField
) const
{
    fValuePtr.set(new dimensionedScalar(dimless, scalar(1)));
}


void Foam::surfaceFunctions::surface::dynamicDiffusivityCoeffs
(
    const LagrangianSubMesh& subMesh,
    autoPtr<dimensionedScalar>& fValuePtr,
    tmp<LagrangianSubScalarField>& tfField
) const
{
    fValuePtr.set(new dimensionedScalar(dimless, scalar(1)));
}


void Foam::surfaceFunctions::surface::kinematicDiffusivityCoeffs
(
    const LagrangianSubMesh& subMesh,
    autoPtr<dimensionedScalar>& fValuePtr,
    tmp<LagrangianSubScalarField>& tfField
) const
{
    fValuePtr.set(new dimensionedScalar(dimless, scalar(1)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFunctions::surface::surface(const LagrangianMesh& mesh)
:
    surfaceFunction(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFunctions::surface::~surface()
{}


// ************************************************************************* //
