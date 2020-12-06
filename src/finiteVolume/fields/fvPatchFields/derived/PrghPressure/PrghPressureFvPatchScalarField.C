/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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

#include "PrghPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PressureFvPatchScalarField>
Foam::PrghPressureFvPatchScalarField<PressureFvPatchScalarField>::
PrghPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    PressureFvPatchScalarField(p, iF),
    rhoName_("rho")
{}


template<class PressureFvPatchScalarField>
Foam::PrghPressureFvPatchScalarField<PressureFvPatchScalarField>::
PrghPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    PressureFvPatchScalarField(p, iF, dict),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{}


template<class PressureFvPatchScalarField>
Foam::PrghPressureFvPatchScalarField<PressureFvPatchScalarField>::
PrghPressureFvPatchScalarField
(
    const PrghPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    PressureFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_)
{}


template<class PressureFvPatchScalarField>
Foam::PrghPressureFvPatchScalarField<PressureFvPatchScalarField>::
PrghPressureFvPatchScalarField
(
    const PrghPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    PressureFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PressureFvPatchScalarField>
void Foam::PrghPressureFvPatchScalarField<PressureFvPatchScalarField>::
updateCoeffs()
{
    if (PressureFvPatchScalarField::updated())
    {
        return;
    }

    PressureFvPatchScalarField::updateCoeffs();

    const scalarField& rhop = this->patch().template
        lookupPatchField<volScalarField, scalar>
        (
            rhoName_
        );

    const uniformDimensionedVectorField& g =
        this->db().template lookupObject<uniformDimensionedVectorField>("g");

    const uniformDimensionedScalarField& hRef =
        this->db().template lookupObject<uniformDimensionedScalarField>("hRef");

    const dimensionedScalar ghRef(- mag(g)*hRef);

    this->operator==
    (
        *this - rhop*((g.value() & this->patch().Cf()) - ghRef.value())
    );
}


template<class PressureFvPatchScalarField>
void Foam::PrghPressureFvPatchScalarField<PressureFvPatchScalarField>::
write(Ostream& os) const
{
    PressureFvPatchScalarField::write(os);
    writeEntry(os, "rho", rhoName_);
}


// ************************************************************************* //
