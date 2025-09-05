/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "meanInletOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meanInletOutletFvPatchScalarField::meanInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF)
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", iF.dimensions(), dict, p.size())
        );
        this->refValue() = *this;
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::meanInletOutletFvPatchScalarField::meanInletOutletFvPatchScalarField
(
    const meanInletOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::meanInletOutletFvPatchScalarField::meanInletOutletFvPatchScalarField
(
    const meanInletOutletFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meanInletOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar averagePsi
    (
        min
        (
            max
            (
                gSum(patch().magSf()*patchInternalField())
               /gSum(patch().magSf()),
                0
            ),
            1
        )
    );

    Info << "averagePsi " << averagePsi << endl;
    refValue() = averagePsi;

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::meanInletOutletFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        meanInletOutletFvPatchScalarField
    );
}

// ************************************************************************* //
