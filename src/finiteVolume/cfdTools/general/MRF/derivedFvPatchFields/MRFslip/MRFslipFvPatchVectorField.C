/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "MRFslipFvPatchVectorField.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFslipFvPatchVectorField::MRFslipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    MRFPatchField(dict)
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(Zero);
    }
}


Foam::MRFslipFvPatchVectorField::MRFslipFvPatchVectorField
(
    const MRFslipFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    MRFPatchField(pvf)
{}


Foam::MRFslipFvPatchVectorField::MRFslipFvPatchVectorField
(
    const MRFslipFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    MRFPatchField(pvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFslipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<vectorField> nHat(patch().nf());

    // Set the patch velocity to the neighbouring cell value
    vectorField::operator=(patchInternalField());

    // Transform the patch velocity relative to the rotating frame
    makeRelative(*this);

    // Transform the patch velocity which removes the normal component
    vectorField::operator=((*this + transform(I - 2*sqr(nHat), *this))/2);

    // Transform the patch velocity absolute to the rotating frame
    makeAbsolute(*this);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::MRFslipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    MRFPatchField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        MRFslipFvPatchVectorField
    );
}

// ************************************************************************* //
