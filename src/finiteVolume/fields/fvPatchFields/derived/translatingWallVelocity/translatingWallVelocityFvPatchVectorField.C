/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "translatingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(vector::zero)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    U_(ptf.U_)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(dict.lookup("U"))
{
    // Evaluate the wall velocity
    updateCoeffs();
}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& twvpvf
)
:
    fixedValueFvPatchField<vector>(twvpvf),
    U_(twvpvf.U_)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& twvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(twvpvf, iF),
    U_(twvpvf.U_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translatingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Remove the component of U normal to the wall in case the wall is not flat
    const vectorField n(patch().nf());
    vectorField::operator=(U_ - n*(n & U_));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::translatingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("U") << U_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        translatingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
