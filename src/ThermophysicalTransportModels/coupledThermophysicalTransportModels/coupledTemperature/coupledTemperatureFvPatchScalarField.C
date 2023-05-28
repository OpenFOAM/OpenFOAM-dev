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

#include "coupledTemperatureFvPatchScalarField.H"
#include "thermophysicalTransportModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledTemperatureFvPatchScalarField::getThis
(
    tmp<scalarField>& kappa,
    tmp<scalarField>& sumKappaTByDelta,
    tmp<scalarField>& sumKappaByDeltaNbr,
    scalarField& sumq,
    tmp<scalarField>& qByKappa
) const
{
    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    kappa = ttm.kappaEff(patch().index());

    qByKappa = sumq/kappa();

    sumq = 0;

    tmp<scalarField> qCorr(ttm.qCorr(patch().index()));

    if (qCorr.valid())
    {
        sumq += qCorr;
    }
}


void Foam::coupledTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& sumKappaTByDeltaNbr,
    tmp<scalarField>& sumKappaByDeltaNbr,
    tmp<scalarField>& qNbr
) const
{
    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    sumKappaByDeltaNbr = ttm.kappaEff(patch().index())*patch().deltaCoeffs();

    sumKappaTByDeltaNbr = sumKappaByDeltaNbr()*patchInternalField();

    qNbr = ttm.qCorr(patch().index());
}


void Foam::coupledTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& TrefNbr,
    tmp<scalarField>& qNbr
) const
{
    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>
        (
            internalField().name()
        );

    TrefNbr = Tp;

    qNbr = ttm.qCorr(patch().index());
}


void Foam::coupledTemperatureFvPatchScalarField::add
(
    tmp<scalarField>& result,
    const tmp<scalarField>& field
) const
{
    if (result.valid())
    {
        result.ref() += field;
    }
    else
    {
        if (field.isTmp())
        {
            result = field;
        }
        else
        {
            result = field().clone();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    qs_(),
    Qs_(0),
    wallKappaByDelta_(0)
{
    mappedPatchBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBase::from::differentPatch
    );

    if (dict.found("thicknessLayers"))
    {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, i)
            {
                wallKappaByDelta_ += thicknessLayers_[i]/kappaLayers_[i];
            }
            wallKappaByDelta_ = 1/wallKappaByDelta_;
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

        qs_ = new scalarField("qs", dict, p.size());
    }
    else if (dict.found("Qs"))
    {
        Qs_ = dict.lookup<scalar>("Qs");
        qs_ = new scalarField(p.size(), Qs_/gSum(patch().magSf()));
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
        refGrad() = 0;
        valueFraction() = 1;
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
    qs_
    (
        psf.qs_.valid()
      ? mapper(psf.qs_()).ptr()
      : nullptr
    ),
    Qs_(psf.Qs_),
    wallKappaByDelta_(psf.wallKappaByDelta_)
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
    qs_(psf.qs_, false),
    Qs_(psf.Qs_),
    wallKappaByDelta_(psf.wallKappaByDelta_)
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

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr =
        refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    const fvPatchScalarField& TpNbr =
        patchNbr.lookupPatchField<volScalarField, scalar>(TnbrName_);

    if (!isA<coupledTemperatureFvPatchScalarField>(TpNbr))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << this->patch().name() << " is of type "
            << coupledTemperatureFvPatchScalarField::typeName
            << endl << "The neighbouring patch field "
            << internalField().name() << " on "
            << patchNbr.name() << " is required to be the same, but is "
            << "currently of type " << TpNbr.type() << exit(FatalError);
    }

    const coupledTemperatureFvPatchScalarField& coupledTemperatureNbr =
        refCast<const coupledTemperatureFvPatchScalarField>(TpNbr);

    scalarField sumq(size(), 0);

    if (qs_.valid())
    {
        sumq += qs_();
    }

    if (qrName_ != "none")
    {
        sumq += patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    if (qrNbrName_ != "none")
    {
        sumq += mpp.fromNeighbour
        (
            patchNbr.lookupPatchField<volScalarField, scalar>(qrNbrName_)
        );
    }

    tmp<scalarField> kappa;
    tmp<scalarField> sumKappaTByDelta;
    tmp<scalarField> sumKappaByDelta;
    tmp<scalarField> qByKappa;

    // q = alpha.this*sumq
    getThis(kappa, sumKappaTByDelta, sumKappaByDelta, sumq, qByKappa);

    // Add neighbour contributions
    {
        tmp<scalarField> sumKappaTByDeltaNbr;
        tmp<scalarField> sumKappaByDeltaNbr;
        tmp<scalarField> qNbr;

        if (wallKappaByDelta_ == 0)
        {
            coupledTemperatureNbr.getNbr
            (
                sumKappaTByDeltaNbr,
                sumKappaByDeltaNbr,
                qNbr
            );

            add(sumKappaTByDelta, mpp.fromNeighbour(sumKappaTByDeltaNbr));
            add(sumKappaByDelta, mpp.fromNeighbour(sumKappaByDeltaNbr));
        }
        else
        {
            // Get the neighbour wall temperature and flux correction
            tmp<scalarField> TwNbr;
            coupledTemperatureNbr.getNbr(TwNbr, qNbr);

            add(sumKappaByDelta, scalarField(size(), wallKappaByDelta_));
            add(sumKappaTByDelta, wallKappaByDelta_*mpp.fromNeighbour(TwNbr));
        }

        if (qNbr.valid())
        {
            sumq += mpp.fromNeighbour(qNbr);
        }
    }

    this->valueFraction() =
        sumKappaByDelta()/(kappa()*patch().deltaCoeffs() + sumKappaByDelta());

    this->refValue() = (sumKappaTByDelta() + sumq)/sumKappaByDelta();

    this->refGrad() = qByKappa;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa()*patch().magSf()*snGrad());

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

    writeEntryIfDifferent<word>(os, "Tnbr", "T", TnbrName_);
    writeEntryIfDifferent<word>(os, "qrNbr", "none", qrNbrName_);
    writeEntryIfDifferent<word>(os, "qr", "none", qrName_);

    if (Qs_ != 0)
    {
        writeEntry(os, "Qs", Qs_);
    }
    else if (qs_.valid())
    {
        writeEntry(os, "qs", qs_());
    }

    if (thicknessLayers_.size())
    {
        writeEntry(os, "thicknessLayers", thicknessLayers_);
        writeEntry(os, "kappaLayers", kappaLayers_);
    }
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
