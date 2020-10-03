/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "filmHeightInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    deltafName_("delta")
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    deltafName_(ptf.deltafName_)
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    deltafName_(dict.lookupOrDefault<word>("delta", "delta"))
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fhivpvf, iF),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const fvPatchField<scalar>& deltafp =
        patch().lookupPatchField<volScalarField, scalar>(deltafName_);

    vectorField n(patch().nf());
    const scalarField& magSf = patch().magSf();

    operator==(n*phip/(rhop*magSf*deltafp + rootVSmall));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmHeightInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "delta", "delta", deltafName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        filmHeightInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
