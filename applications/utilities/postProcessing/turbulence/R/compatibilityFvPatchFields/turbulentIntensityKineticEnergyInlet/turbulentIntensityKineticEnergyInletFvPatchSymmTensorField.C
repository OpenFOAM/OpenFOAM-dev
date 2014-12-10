/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "turbulentIntensityKineticEnergyInletFvPatchSymmTensorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentIntensityKineticEnergyInletFvPatchSymmTensorField::
turbulentIntensityKineticEnergyInletFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(p, iF)
{}


Foam::turbulentIntensityKineticEnergyInletFvPatchSymmTensorField::
turbulentIntensityKineticEnergyInletFvPatchSymmTensorField
(
    const turbulentIntensityKineticEnergyInletFvPatchSymmTensorField& ptf,
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchSymmTensorField(ptf, p, iF, mapper)
{}


Foam::turbulentIntensityKineticEnergyInletFvPatchSymmTensorField::
turbulentIntensityKineticEnergyInletFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchSymmTensorField(p, iF)
{

    fvPatchSymmTensorField::operator=(symmTensorField("value", dict, p.size()));
}


Foam::turbulentIntensityKineticEnergyInletFvPatchSymmTensorField::
turbulentIntensityKineticEnergyInletFvPatchSymmTensorField
(
    const turbulentIntensityKineticEnergyInletFvPatchSymmTensorField& ptf
)
:
    fixedValueFvPatchSymmTensorField(ptf)
{}


Foam::turbulentIntensityKineticEnergyInletFvPatchSymmTensorField::
turbulentIntensityKineticEnergyInletFvPatchSymmTensorField
(
    const turbulentIntensityKineticEnergyInletFvPatchSymmTensorField& ptf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentIntensityKineticEnergyInletFvPatchSymmTensorField::write
(
    Ostream& os
) const
{
    fvPatchSymmTensorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchSymmTensorField,
        turbulentIntensityKineticEnergyInletFvPatchSymmTensorField
    );
}

// ************************************************************************* //
