/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "surfaceNormalUniformFixedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceNormalUniformFixedValueFvPatchVectorField::
surfaceNormalUniformFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    uniformValue_(Function1<scalar>::New("uniformValue", dict))
{
    this->evaluate();
}


Foam::surfaceNormalUniformFixedValueFvPatchVectorField::
surfaceNormalUniformFixedValueFvPatchVectorField
(
    const surfaceNormalUniformFixedValueFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper, false), // Don't map
    uniformValue_(pvf.uniformValue_, false)
{
    // Evaluate since value not mapped
    this->evaluate();
}


Foam::surfaceNormalUniformFixedValueFvPatchVectorField::
surfaceNormalUniformFixedValueFvPatchVectorField
(
    const surfaceNormalUniformFixedValueFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    uniformValue_(pvf.uniformValue_, false)
{
    // Evaluate the profile if defined
    if (pvf.uniformValue_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceNormalUniformFixedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().userTimeValue();
    fvPatchVectorField::operator=(uniformValue_->value(t)*patch().nf());

    fvPatchVectorField::updateCoeffs();
}


void Foam::surfaceNormalUniformFixedValueFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, uniformValue_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        surfaceNormalUniformFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //
