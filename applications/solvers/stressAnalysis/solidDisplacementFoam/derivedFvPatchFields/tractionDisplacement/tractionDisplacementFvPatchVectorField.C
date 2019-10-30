/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "tractionDisplacementFvPatchVectorField.H"
#include "solidDisplacementThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
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


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(mapper(tdpvf.traction_)),
    pressure_(mapper(tdpvf.pressure_))
{}


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
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


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


Foam::tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    m(traction_, traction_);
    m(pressure_, pressure_);
}


void Foam::tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


void Foam::tractionDisplacementFvPatchVectorField::updateCoeffs()
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
    const scalarField threeK
    (
        thermo.planeStress()
      ? E/(1 - nu)
      : E/(1 - 2*nu)
    );

    const scalarField twoMuLambda(2*mu + lambda);

    const vectorField n(patch().nf());

    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    gradient() =
    (
        (traction_ - pressure_*n)
      + twoMuLambda*fvPatchField<vector>::snGrad() - (n & sigmaD)
    )/twoMuLambda;

    if (thermo.thermalStress())
    {
        const scalarField& alphav = thermo.alphav(patchi);

        gradient() +=
            n*threeK*alphav*thermo.T().boundaryField()[patchi]/twoMuLambda;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "traction", traction_);
    writeEntry(os, "pressure", pressure_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tractionDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
