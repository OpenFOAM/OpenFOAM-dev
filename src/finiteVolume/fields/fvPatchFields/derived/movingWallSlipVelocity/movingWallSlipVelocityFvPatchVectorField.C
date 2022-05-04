/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "movingWallSlipVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF)
{
    const vectorField n(p.nf());

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = sqr(n);
}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    const vectorField n(p.nf());

    refValue() = sqr(n) & *this;
    refGrad() = Zero;
    valueFraction() = sqr(n);
}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    if (mesh.moving())
    {
        const fvPatch& p = patch();

        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());

        const scalarField phip(fvc::meshPhi(U, p.index()));

        const vectorField n(p.nf());
        const scalarField& magSf = p.magSf();

        tmp<scalarField> Un = phip/(magSf + vSmall);

        refValue() = n*Un;
        refGrad() = Zero;
        valueFraction() = sqr(n);
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::movingWallSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingWallSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
