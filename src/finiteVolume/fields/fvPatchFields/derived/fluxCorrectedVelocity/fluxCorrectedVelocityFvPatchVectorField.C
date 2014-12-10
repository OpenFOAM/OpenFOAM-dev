/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "fluxCorrectedVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxCorrectedVelocityFvPatchVectorField::
fluxCorrectedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::fluxCorrectedVelocityFvPatchVectorField::
fluxCorrectedVelocityFvPatchVectorField
(
    const fluxCorrectedVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::fluxCorrectedVelocityFvPatchVectorField::
fluxCorrectedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    fvPatchVectorField::operator=(patchInternalField());
}


Foam::fluxCorrectedVelocityFvPatchVectorField::
fluxCorrectedVelocityFvPatchVectorField
(
    const fluxCorrectedVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(fcvpvf, iF),
    phiName_(fcvpvf.phiName_),
    rhoName_(fcvpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxCorrectedVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    zeroGradientFvPatchVectorField::evaluate();

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    const vectorField n(patch().nf());
    const Field<scalar>& magS = patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(*this - n*(n & *this) + n*phip/magS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(*this - n*(n & *this) + n*phip/(rhop*magS));
    }
    else
    {
        FatalErrorIn
        (
            "fluxCorrectedVelocityFvPatchVectorField::evaluate()"
        )
            << "dimensions of phi are incorrect\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


void Foam::fluxCorrectedVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fluxCorrectedVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
