/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "patchPatchDist.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::flowRateInletVelocityFvPatchVectorField::setWallDist()
{
    if (profile_.valid())
    {
        const labelHashSet otherPatchIDs
        (
            patch().patch().boundaryMesh().findPatchIDs<wallPolyPatch>()
        );

        const patchPatchDist pwd(patch().patch(), otherPatchIDs);

        y_ = pwd/gMax(pwd);
    }

    area_ = gSum(patch().magSf());
}


Foam::tmp<Foam::scalarField>
Foam::flowRateInletVelocityFvPatchVectorField::profile()
{
    if (profile_.valid())
    {
        return profile_->value(y_);
    }
    else
    {
        return tmp<scalarField>(new scalarField(size(), scalar(1)));
    }
}


template<class ScaleType, class AlphaType, class RhoType>
void Foam::flowRateInletVelocityFvPatchVectorField::updateValues
(
    const ScaleType& scale,
    const AlphaType& alpha,
    const RhoType& rho
)
{
    const scalarField profile(this->profile());

    const scalar avgU =
      -(scale*flowRate_->value(db().time().userTimeValue()))
       /gSum(alpha*rho*profile*patch().magSf());

    operator==(avgU*profile*patch().nf());
}


template<class AlphaType>
void Foam::flowRateInletVelocityFvPatchVectorField::updateValues
(
    const AlphaType& alpha
)
{
    if (meanVelocity_)
    {
        updateValues(area_, alpha, one());
    }
    else if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one(), alpha, one());
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(one(), alpha, rhop);
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

            updateValues(one(), alpha, rhoInlet_);
        }
    }
}


bool Foam::flowRateInletVelocityFvPatchVectorField::canEvaluate()
{
    return
        Pstream::parRun()
     || !patch().boundaryMesh().mesh().time().processorCase();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    rhoName_("rho"),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -vGreat)),
    alphaName_(dict.lookupOrDefault<word>("alpha", word::null))
{
    if (dict.found("meanVelocity"))
    {
        meanVelocity_ = true;
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("meanVelocity", dict);
    }
    else if (dict.found("volumetricFlowRate"))
    {
        meanVelocity_ = false;
        volumetric_ = true;
        flowRate_ = Function1<scalar>::New("volumetricFlowRate", dict);
    }
    else if (dict.found("massFlowRate"))
    {
        meanVelocity_ = false;
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("massFlowRate", dict);
        rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Please supply 'meanVelocity', 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
    }

    if (dict.found("profile"))
    {
        profile_ = Function1<scalar>::New("profile", dict);
    }

    if (canEvaluate())
    {
        setWallDist();
    }

    if (!canEvaluate() || dict.found("value"))
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
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_, false),
    profile_(ptf.profile_, false),
    meanVelocity_(ptf.meanVelocity_),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    alphaName_(ptf.alphaName_),
    y_
    (
        profile_.valid() && canEvaluate()
      ? mapper(ptf.y_)
      : tmp<scalarField>(new scalarField())
    ),
    area_(NaN)
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
    profile_(ptf.profile_, false),
    meanVelocity_(ptf.meanVelocity_),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    alphaName_(ptf.alphaName_),
    y_(ptf.y_),
    area_(ptf.area_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    const flowRateInletVelocityFvPatchVectorField& tiptf =
        refCast<const flowRateInletVelocityFvPatchVectorField>(ptf);

    if (profile_.valid() && canEvaluate())
    {
        mapper(y_, tiptf.y_);
    }
}


void Foam::flowRateInletVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    const flowRateInletVelocityFvPatchVectorField& tiptf =
        refCast<const flowRateInletVelocityFvPatchVectorField>(ptf);

    if (profile_.valid() && canEvaluate())
    {
        y_.reset(tiptf.y_);
    }
}


void Foam::flowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!canEvaluate())
    {
        FatalErrorInFunction
            << "Cannot evaluate flow rate on a non-parallel processor case"
            << exit(FatalError);
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
    if (profile_.valid())
    {
        writeEntry(os, profile_());
    }
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoInlet", -vGreat, rhoInlet_);
    }
    writeEntryIfDifferent<word>(os, "alpha", word::null, alphaName_);
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
