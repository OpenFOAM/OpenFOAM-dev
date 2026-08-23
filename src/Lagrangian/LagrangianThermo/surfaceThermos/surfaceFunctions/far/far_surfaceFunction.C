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

#include "far_surfaceFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFunctions
{
    defineTypeNameAndDebug(far, 0);
    addToRunTimeSelectionTable(surfaceFunction, far, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceFunctions::far::thermalConductivityCoeffs
(
    const LagrangianSubMesh& subMesh,
    autoPtr<dimensionedScalar>& fValuePtr,
    tmp<LagrangianSubScalarField>& tfField
) const
{
    fValuePtr.set(new dimensionedScalar(dimless, scalar(0)));
}


void Foam::surfaceFunctions::far::dynamicDiffusivityCoeffs
(
    const LagrangianSubMesh& subMesh,
    autoPtr<dimensionedScalar>& fValuePtr,
    tmp<LagrangianSubScalarField>& tfField
) const
{
    fValuePtr.set(new dimensionedScalar(dimless, scalar(0)));
}


void Foam::surfaceFunctions::far::kinematicDiffusivityCoeffs
(
    const LagrangianSubMesh& subMesh,
    autoPtr<dimensionedScalar>& fValuePtr,
    tmp<LagrangianSubScalarField>& tfField
) const
{
    fValuePtr.set(new dimensionedScalar(dimless, scalar(0)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFunctions::far::far(const LagrangianMesh& mesh)
:
    surfaceFunction(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFunctions::far::~far()
{}


// ************************************************************************* //
