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

#include "flowRateInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

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
    rhoInlet_(0.0),
    alphaName_(word::null),
    extrapolateProfile_(false)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -vGreat)),
    alphaName_(dict.lookupOrDefault<word>("alpha", word::null)),
    extrapolateProfile_
    (
        dict.lookupOrDefault<Switch>("extrapolateProfile", false)
    )
{
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ = Function1<scalar>::New("volumetricFlowRate", dict);
        rhoName_ = "rho";
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("massFlowRate", dict);
        rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
        FatalIOErrorInFunction
        (
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
        evaluate(Pstream::commsTypes::blocking);
    }
}


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
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    alphaName_(ptf.alphaName_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    alphaName_(ptf.alphaName_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class AlphaType, class RhoType>
void Foam::flowRateInletVelocityFvPatchVectorField::updateValues
(
    const AlphaType& alpha,
    const RhoType& rho
)
{
    const scalar t = db().time().timeOutputValue();

    const vectorField n(patch().nf());

    if (extrapolateProfile_)
    {
        vectorField Up(this->patchInternalField());

        // Patch normal extrapolated velocity
        scalarField nUp(n & Up);

        // Remove the normal component of the extrapolate patch velocity
        Up -= nUp*n;

        // Remove any reverse flow
        nUp = min(nUp, scalar(0));

        const scalar flowRate = flowRate_->value(t);
        const scalar estimatedFlowRate =
            -gSum(alpha*rho*(this->patch().magSf()*nUp));

        if (estimatedFlowRate/flowRate > 0.5)
        {
            nUp *= mag(flowRate)/mag(estimatedFlowRate);
        }
        else
        {
            nUp -=
                (flowRate - estimatedFlowRate)
               /gSum(alpha*rho*patch().magSf());
        }

        // Add the corrected normal component of velocity to the patch velocity
        Up += nUp*n;

        // Correct the patch velocity
        this->operator==(Up);
    }
    else
    {
        const scalar avgU =
            -flowRate_->value(t)/gSum(alpha*rho*patch().magSf());
        operator==(avgU*n);
    }
}


template<class AlphaType>
void Foam::flowRateInletVelocityFvPatchVectorField::updateValues
(
    const AlphaType& alpha
)
{
    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(alpha, one());
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(alpha, rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(alpha, rhoInlet_);
        }
    }
}


void Foam::flowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (alphaName_ != word::null)
    {
        const fvPatchField<scalar>& alphap =
            patch().lookupPatchField<volScalarField, scalar>(alphaName_);

        updateValues(alphap);
    }
    else
    {
        updateValues(one());
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, flowRate_());
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoInlet", -vGreat, rhoInlet_);
    }
    writeEntryIfDifferent<word>(os, "alpha", word::null, alphaName_);
    writeEntry(os, "extrapolateProfile", extrapolateProfile_);
    writeEntry(os, "value", *this);
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
