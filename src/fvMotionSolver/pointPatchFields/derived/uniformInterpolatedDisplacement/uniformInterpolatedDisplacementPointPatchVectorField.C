/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "uniformInterpolatedDisplacementPointPatchVectorField.H"
#include "pointFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformInterpolatedDisplacementPointPatchVectorField::
uniformInterpolatedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    pointInterpolator_
    (
        new dynamicMeshPointInterpolator(iF.mesh().mesh(), dict)
    )
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::uniformInterpolatedDisplacementPointPatchVectorField::
uniformInterpolatedDisplacementPointPatchVectorField
(
    const uniformInterpolatedDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::uniformInterpolatedDisplacementPointPatchVectorField::
uniformInterpolatedDisplacementPointPatchVectorField
(
    const uniformInterpolatedDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformInterpolatedDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Extract back from the internal field
    this->operator==
    (
        patchInternalField(pointInterpolator_->curPointField()())
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::uniformInterpolatedDisplacementPointPatchVectorField::write
(
    Ostream& os
)
const
{
    pointPatchField<vector>::write(os);
    pointInterpolator_->write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        uniformInterpolatedDisplacementPointPatchVectorField
    );
}

// ************************************************************************* //
