/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),
    TnbrName_("undefined-Tnbr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    qs_(p.size()),
    contactRes_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    TnbrName_(dict.lookup("Tnbr")),
    thicknessLayers_(0),
    kappaLayers_(0),
    qs_(p.size(), 0),
    contactRes_(0.0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("thicknessLayers"))
    {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    if (dict.found("qs"))
    {
        if (dict.found("Qs"))
        {
            FatalIOErrorInFunction(dict)
                << "Either qs or Qs should be specified, not both"
                << exit(FatalIOError);
        }

        qs_ = scalarField("qs", dict, p.size());
    }
    else if (dict.found("Qs"))
    {
        qs_ = scalarField
        (
            p.size(),
            dict.lookup<scalar>("Qs")/gSum(patch().magSf())
        );
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    qs_(mapper(ptf.qs_)),
    contactRes_(ptf.contactRes_)
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    temperatureCoupledBase(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    qs_(ptf.qs_),
    contactRes_(ptf.contactRes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    typedef turbulentTemperatureCoupledBaffleMixedFvPatchScalarField thisType;

    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);

    if (!isA<thisType>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << TnbrName_ << " on "
            << nbrPatch.name() << " is required to be the same, but is "
            << "currently of type " << nbrTp.type() << exit(FatalError);
    }

    const scalarField myKDelta(kappa(*this)*patch().deltaCoeffs());

    const thisType& nbrField = refCast<const thisType>(nbrTp);

    const scalarField nbrIntFld
    (
        contactRes_ == 0
      ? mpp.distribute(nbrField.patchInternalField())
      : mpp.distribute(nbrField)
    );
    const scalarField nbrKDelta
    (
        contactRes_ == 0
      ? mpp.distribute(nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs())
      : tmp<scalarField>(new scalarField(contactRes_, nbrField.size()))
    );

    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = qs_/kappa;
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta)

    this->refValue() = nbrIntFld;
    this->refGrad() = qs_/kappa(*this);
    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta);

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "Tnbr", TnbrName_);
    writeEntry(os, "thicknessLayers", thicknessLayers_);
    writeEntry(os, "kappaLayers", kappaLayers_);

    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
