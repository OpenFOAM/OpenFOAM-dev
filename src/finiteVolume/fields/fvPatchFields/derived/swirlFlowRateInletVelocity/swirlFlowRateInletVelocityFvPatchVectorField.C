/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "swirlFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class RhoType>
void Foam::swirlFlowRateInletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const scalar t = db().time().value();
    const scalarField ts(size(), t);

    // Compute geometry
    const vector axisHat = normalised(axis_);
    const vectorField d(patch().Cf() - origin_);
    const vectorField r(d - (axisHat & d)*axisHat);
    const scalarField magR(mag(r));
    const vectorField rHat(normalised(r));

    // Evaluate individual velocity components
    const scalar axialVelocity = flowRate_->value(t)/gSum(rho*patch().magSf());
    const scalarField radialVelocity(radialVelocity_->value(ts, magR));
    tmp<scalarField> tangentialVelocity;
    if (omega_.valid())
    {
        tangentialVelocity = omega_->value(t)*magR;
    }
    else
    {
        tangentialVelocity = tangentialVelocity_->value(ts, magR);
    }

    // Combine components the complete vector velocity
    operator==
    (
        axialVelocity*axisHat
      + radialVelocity*rHat
      + tangentialVelocity*(axisHat ^ rHat)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    origin_
    (
        dict.lookupOrDefault
        (
            "origin",
            dimLength,
            returnReduce(patch().size(), sumOp<label>())
          ? gSum(patch().Cf()*patch().magSf())/gSum(patch().magSf())
          : Zero
        )
    ),
    axis_
    (
        dict.lookupOrDefault
        (
            "axis",
            dimless,
            returnReduce(patch().size(), sumOp<label>())
          ? -gSum(patch().Sf())/gSum(patch().magSf())
          : Zero
        )
    ),
    flowRate_(),
    volumetric_(),
    rhoName_("rho"),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", dimDensity, -vGreat)),
    radialVelocity_
    (
        Function2<scalar>::New
        (
            "radialVelocity",
            db().time().userUnits(),
            dimLength,
            dimVelocity,
            dict
        )
    ),
    omega_(nullptr),
    tangentialVelocity_(nullptr)
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

    if (dict.found("omega") || dict.found("rpm"))
    {
        omega_ = new Function1s::omega(db().time(), dict);
    }
    else if (dict.found("tangentialVelocity"))
    {
        tangentialVelocity_ =
            Function2<scalar>::New
            (
                "tangentialVelocity",
                db().time().userUnits(),
                dimLength,
                dimVelocity,
                dict
            );
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'omega' or 'rpm' or"
            << " 'tangentialVelocity'" << exit(FatalIOError);
    }

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


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    radialVelocity_(ptf.radialVelocity_, false),
    omega_(ptf.omega_, false),
    tangentialVelocity_(ptf.tangentialVelocity_, false)
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    radialVelocity_(ptf.radialVelocity_, false),
    omega_(ptf.omega_, false),
    tangentialVelocity_(ptf.tangentialVelocity_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
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
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::swirlFlowRateInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "axis", axis_);
    writeEntry(os, db().time().userUnits(), unitAny, flowRate_());
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoInlet", -vGreat, rhoInlet_);
    }
    writeEntry
    (
        os,
        db().time().userUnits(),
        dimLength,
        dimVelocity,
        radialVelocity_()
    );
    if (omega_.valid())
    {
        writeEntry(os, omega_());
    }
    else
    {
        writeEntry
        (
            os,
            db().time().userUnits(),
            dimLength,
            dimVelocity,
            tangentialVelocity_()
        );
    }
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       swirlFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
