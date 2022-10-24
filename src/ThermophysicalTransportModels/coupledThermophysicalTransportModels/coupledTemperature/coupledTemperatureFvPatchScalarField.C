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

#include "coupledTemperatureFvPatchScalarField.H"
#include "thermophysicalTransportModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledTemperatureFvPatchScalarField::
coupledTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    qs_(p.size()),
    contactRes_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


Foam::coupledTemperatureFvPatchScalarField::
coupledTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.lookupOrDefault<word>("qrNbr", "none")),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
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


Foam::coupledTemperatureFvPatchScalarField::
coupledTemperatureFvPatchScalarField
(
    const coupledTemperatureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    qs_(mapper(psf.qs_)),
    contactRes_(psf.contactRes_)
{}


Foam::coupledTemperatureFvPatchScalarField::
coupledTemperatureFvPatchScalarField
(
    const coupledTemperatureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    qs_(psf.qs_),
    contactRes_(psf.contactRes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const label patchi = patch().index();

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr =
        refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    typedef coupledTemperatureFvPatchScalarField thisType;

    const fvPatchScalarField& TpNbr =
        patchNbr.lookupPatchField<volScalarField, scalar>(TnbrName_);

    if (!isA<thisType>(TpNbr))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << TnbrName_ << " on "
            << patchNbr.name() << " is required to be the same, but is "
            << "currently of type " << TpNbr.type() << exit(FatalError);
    }

    const thisType& coupledTemperatureNbr = refCast<const thisType>(TpNbr);

    const scalarField TcNbr
    (
        contactRes_ == 0
      ? mpp.distribute(coupledTemperatureNbr.patchInternalField())
      : mpp.distribute(coupledTemperatureNbr)
    );

    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    const thermophysicalTransportModel& ttmNbr =
        patchNbr.boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    const scalarField kappa(ttm.kappaEff(patchi));

    const scalarField KDelta(kappa*patch().deltaCoeffs());

    const scalarField KDeltaNbr
    (
        contactRes_ == 0
      ? mpp.distribute
        (
            ttmNbr.kappaEff(patchiNbr)*patchNbr.deltaCoeffs()
        )
      : tmp<scalarField>(new scalarField(size(), contactRes_))
    );

    scalarField qTot(qs_);

    if (qrName_ != "none")
    {
        qTot += patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    if (qrNbrName_ != "none")
    {
        qTot += mpp.distribute
        (
            patchNbr.lookupPatchField<volScalarField, scalar>(qrNbrName_)
        );
    }

    tmp<scalarField> qCorr(ttm.qCorr(patchi));

    if (qCorr.valid())
    {
        qTot += qCorr;
    }

    tmp<scalarField> qCorrNbr(ttmNbr.qCorr(patchiNbr));

    if (qCorrNbr.valid())
    {
        qTot += mpp.distribute(qCorrNbr);
    }

    // Both sides agree on
    // - temperature : (KDelta*fld + KDeltaNbr*nbrFld)/(KDelta + KDeltaNbr)
    // - gradient    : (temperature - fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedValue and pure fixedGradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = qTot/kappa;
    //    - refValue = neighbour value
    //    - mixFraction = KDeltaNbr / (KDeltaNbr + KDelta)

    this->valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
    this->refValue() = TcNbr;
    this->refGrad() = qTot/kappa;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.nbrMesh().name() << ':'
            << patchNbr.name() << ':'
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


void Foam::coupledTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "Tnbr", TnbrName_);
    writeEntry(os, "qrNbr", qrNbrName_);
    writeEntry(os, "qr", qrName_);
    writeEntry(os, "thicknessLayers", thicknessLayers_);
    writeEntry(os, "kappaLayers", kappaLayers_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField
    );


    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField,
        patch,
        turbulentTemperatureCoupledBaffleMixed,
        "compressible::turbulentTemperatureCoupledBaffleMixed"
    );

    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField,
        patchMapper,
        turbulentTemperatureCoupledBaffleMixed,
        "compressible::turbulentTemperatureCoupledBaffleMixed"
    );

    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField,
        dictionary,
        turbulentTemperatureCoupledBaffleMixed,
        "compressible::turbulentTemperatureCoupledBaffleMixed"
    );


    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField,
        patch,
        turbulentTemperatureRadCoupledMixed,
        "compressible::turbulentTemperatureRadCoupledMixed"
    );

    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField,
        patchMapper,
        turbulentTemperatureRadCoupledMixed,
        "compressible::turbulentTemperatureRadCoupledMixed"
    );

    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        coupledTemperatureFvPatchScalarField,
        dictionary,
        turbulentTemperatureRadCoupledMixed,
        "compressible::turbulentTemperatureRadCoupledMixed"
    );
}


// ************************************************************************* //
