/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "outletPhaseMeanVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletPhaseMeanVelocityFvPatchVectorField::
outletPhaseMeanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF, dict, false),
    UnMean_(Function1<scalar>::New("UnMean", dict)),
    alphaName_(dict.lookup("alpha"))
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


Foam::outletPhaseMeanVelocityFvPatchVectorField::
outletPhaseMeanVelocityFvPatchVectorField
(
    const outletPhaseMeanVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    UnMean_(ptf.UnMean_, false),
    alphaName_(ptf.alphaName_)
{}


Foam::outletPhaseMeanVelocityFvPatchVectorField::
outletPhaseMeanVelocityFvPatchVectorField
(
    const outletPhaseMeanVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    UnMean_(ptf.UnMean_, false),
    alphaName_(ptf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletPhaseMeanVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().userTimeValue();

    scalarField alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);

    alphap = max(alphap, scalar(0));
    alphap = min(alphap, scalar(1));

    // Get the patchInternalField (zero-gradient field)
    vectorField Uzg(patchInternalField());

    // Calculate the phase mean zero-gradient normal velocity
    const scalar UnzgMean =
        gSum(alphap*(patch().Sf() & Uzg))
       /gSum(alphap*patch().magSf());

    // Set the refValue and valueFraction to adjust the boundary field
    // such that the phase mean is UnMean_
    const scalar UnMean = UnMean_->value(t);
    if (UnzgMean >= UnMean)
    {
        refValue() = Zero;
        valueFraction() = 1.0 - UnMean/UnzgMean;
    }
    else
    {
        refValue() = (UnMean + UnzgMean)*patch().nf();
        valueFraction() = 1.0 - UnzgMean/UnMean;
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::outletPhaseMeanVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    writeEntry(os, UnMean_());
    writeEntry(os, "alpha", alphaName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       outletPhaseMeanVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
