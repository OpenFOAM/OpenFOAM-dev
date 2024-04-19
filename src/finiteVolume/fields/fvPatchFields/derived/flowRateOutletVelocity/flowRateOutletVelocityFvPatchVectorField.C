/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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

#include "flowRateOutletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class RhoType>
void Foam::flowRateOutletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const vectorField n(patch().nf());

    // Extrapolate patch velocity
    vectorField Up(this->patchInternalField());

    // Patch normal extrapolated velocity
    scalarField nUp(n & Up);

    // Remove the normal component of the extrapolate patch velocity
    Up -= nUp*n;

    // Remove any reverse flow
    nUp = max(nUp, scalar(0));

    const scalar flowRate = flowRate_->value(db().time().value());
    const scalar estimatedFlowRate = gSum(rho*(this->patch().magSf()*nUp));

    if (estimatedFlowRate/flowRate > 0.5)
    {
        nUp *= (mag(flowRate)/mag(estimatedFlowRate));
    }
    else
    {
        nUp += ((flowRate - estimatedFlowRate)/gSum(rho*patch().magSf()));
    }

    // Add the corrected normal component of velocity to the patch velocity
    Up += nUp*n;

    // Correct the patch velocity
    this->operator==(Up);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    flowRate_(),
    volumetric_(),
    rhoName_("rho"),
    rhoOutlet_(dict.lookupOrDefault<scalar>("rhoOutlet", dimDensity, -vGreat))
{
    if (dict.found("volumetricFlowRate"))
    {
        flowRate_ =
            Function1<scalar>::New
            (
                "volumetricFlowRate",
                db().time().userUnits(),
                dimVolumetricFlux,
                dict
            );
        volumetric_ = true;
    }
    else if (dict.found("massFlowRate"))
    {
        flowRate_ =
            Function1<scalar>::New
            (
                "massFlowRate",
                db().time().userUnits(),
                dimMassFlux,
                dict
            );
        volumetric_ = false;
        rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
    }

    // Value field required if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const flowRateOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_)
{}


Foam::flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const flowRateOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one());
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoOutlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoOutlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoOutlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, db().time().userUnits(), unitAny, flowRate_());
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoOutlet", -vGreat, rhoOutlet_);
    }
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateOutletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
