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

#include "tractionDisplacementCorrectionFvPatchVectorField.H"
#include "solidDisplacementThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), Zero),
    pressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const tractionDisplacementCorrectionFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(mapper(tdpvf.traction_)),
    pressure_(mapper(tdpvf.pressure_))
{}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const tractionDisplacementCorrectionFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementCorrectionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    m(traction_, traction_);
    m(pressure_, pressure_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementCorrectionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementCorrectionFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementCorrectionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementCorrectionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const solidDisplacementThermo& thermo =
        db().lookupObject<solidDisplacementThermo>
        (
            solidDisplacementThermo::dictName
        );

    const scalarField& E = thermo.E(patchi);
    const scalarField& nu = thermo.nu(patchi);

    const scalarField mu(E/(2.0*(1.0 + nu)));
    const scalarField lambda
    (
        thermo.planeStress()
      ? nu*E/((1 + nu)*(1 - nu))
      : nu*E/((1 + nu)*(1 - 2*nu))
    );

    const vectorField n(patch().nf());

    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    const fvPatchField<tensor>& sigmaExp =
        patch().lookupPatchField<volTensorField, tensor>("sigmaExp");

    gradient() =
    (
        (traction_ + pressure_*n) - (n & (sigmaD + sigmaExp))
    )/(2.0*mu + lambda);

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void tractionDisplacementCorrectionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "traction", traction_);
    writeEntry(os, "pressure", pressure_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementCorrectionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
