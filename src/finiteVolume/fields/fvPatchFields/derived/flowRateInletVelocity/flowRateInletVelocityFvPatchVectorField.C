/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "flowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
    volumetric_(false),
    rhoName_("rho"),
    rhoInlet_(0.0)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_().clone().ptr()),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -VGREAT))
{
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ = DataEntry<scalar>::New("volumetricFlowRate", dict);
        rhoName_ = "rho";
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = DataEntry<scalar>::New("massFlowRate", dict);
        rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
        FatalIOErrorIn
        (
            "flowRateInletVelocityFvPatchVectorField::"
            "flowRateInletVelocityFvPatchVectorField"
            "(const fvPatch&, const DimensionedField<vector, volMesh>&,"
            " const dictionary&)",
            dict
        )   << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
    }

    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::blocking);
    }
}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_().clone().ptr()),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_().clone().ptr()),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    // a simpler way of doing this would be nice
    const scalar avgU = -flowRate_->value(t)/gSum(patch().magSf());

    tmp<vectorField> n = patch().nf();

    if (volumetric_ || rhoName_ == "none")
    {
        // volumetric flow-rate or density not given
        operator==(n*avgU);
    }
    else
    {
        // mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            operator==(n*avgU/rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorIn
                (
                    "flowRateInletVelocityFvPatchVectorField::updateCoeffs()"
                )   << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }
            operator==(n*avgU/rhoInlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoInlet", -VGREAT, rhoInlet_);
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
