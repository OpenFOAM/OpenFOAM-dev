/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
#include "FunctionalDimensionedField.H"
#include "mappedFvPatchBaseBase.H"
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
        patch().mesh()
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
        patch().mesh()
       .lookupType<thermophysicalTransportModel>();

    sumKappaByDeltaNbr = ttm.kappaEff(patch().index())*patch().deltaCoeffs();

    sumKappaTByDeltaNbr = sumKappaByDeltaNbr()*patchInternalField();

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
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.lookupOrDefault<word>("qrNbr", "none")),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    qrRelax_(dict.lookupOrDefault<scalar>("qrRelaxation", units::fraction, 1)),
    qrPrevious_
    (
        qrName_ != word::null
      ? dict.found("qrPrevious")
      ? scalarField
        (
            "qrPrevious",
            dimensions::heatFluxDensity,
            dict,
            p.size()
        )
      : scalarField(p.size(), 0)
      : scalarField()
    ),
    h_
    (
        dict.found("h")
      ? new FunctionalDimensionedField<scalar, fvPatch>
        (
            iF.name(),
            "h",
            p,
            dimensions::heatFluxDensity/dimensions::temperature,
            dict
        )
      : nullptr
    ),
    qs_(),
    Qs_(0)
{
    mappedPatchBaseBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBaseBase::from::differentPatch
    );

    if (dict.found("qs"))
    {
        if (dict.found("Qs"))
        {
            FatalIOErrorInFunction(dict)
                << "Either qs or Qs should be specified, not both"
                << exit(FatalIOError);
        }

        qs_ = new scalarField
        (
            "qs",
            dimensions::power/dimensions::time,
            dict,
            p.size()
       );
    }
    else if (dict.found("Qs"))
    {
        Qs_ = dict.lookup<scalar>("Qs");
        qs_ = new scalarField(p.size(), Qs_/gSum(patch().magSf()));
    }

    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", iF.dimensions(), dict, p.size());
        refGrad() =
            scalarField
            (
                "refGradient",
                iF.dimensions()/dimensions::length,
                dict,
                p.size()
            );
        valueFraction() =
            scalarField("valueFraction", units::fraction, dict, p.size());
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
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    qrRelax_(psf.qrRelax_),
    qrPrevious_
    (
        qrName_ != word::null
      ? mapper(psf.qrPrevious_)()
      : scalarField()
    ),
    h_
    (
        psf.h_.valid()
      ? new FunctionalDimensionedField<scalar, fvPatch>
        (
            psf.h_(),
            p,
            mapper
        )
      : nullptr
    ),
    qs_(psf.qs_.valid() ? new scalarField(p.size()) : nullptr),
    Qs_(psf.Qs_)
{
    map(psf, mapper);
}


Foam::coupledTemperatureFvPatchScalarField::
coupledTemperatureFvPatchScalarField
(
    const coupledTemperatureFvPatchScalarField& psf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    qrRelax_(psf.qrRelax_),
    qrPrevious_(psf.qrPrevious_),
    h_
    (
        psf.h_.valid()
      ? new FunctionalDimensionedField<scalar, fvPatch>
        (
            psf.h_()
        )
      : nullptr
    ),
    qs_(psf.qs_, false),
    Qs_(psf.Qs_)
{}


Foam::tmp<Foam::fvPatchScalarField>
Foam::coupledTemperatureFvPatchScalarField::clone
(
    const DimensionedField<scalar, fvMesh>& iF
) const
{
    return tmp<fvPatchScalarField>
    (
        new coupledTemperatureFvPatchScalarField(*this, iF)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledTemperatureFvPatchScalarField::
~coupledTemperatureFvPatchScalarField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledTemperatureFvPatchScalarField::map
(
    const coupledTemperatureFvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    // Unmapped faces are considered zero-gradient/adiabatic
    mapper(*this, ptf, [&](){ return patchInternalField(); });
    mapper(refValue(), ptf.refValue(), [&](){ return patchInternalField(); });
    mapper(refGrad(), ptf.refGrad(), scalar(0));
    mapper(valueFraction(), ptf.valueFraction(), scalar(0));

    // Map the heat flux, if present
    if (ptf.qs_.valid())
    {
        mapper(qs_(), ptf.qs_());
    }
}


void Foam::coupledTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    const coupledTemperatureFvPatchScalarField& ctptf =
        refCast<const coupledTemperatureFvPatchScalarField>(ptf);

    map(ctptf, mapper);

    if (h_.valid())
    {
        h_->map(ctptf.h_(), mapper);
    }
}


void Foam::coupledTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const coupledTemperatureFvPatchScalarField& tiptf =
        refCast<const coupledTemperatureFvPatchScalarField>(ptf);

    if (tiptf.qs_.valid())
    {
        qs_().reset(tiptf.qs_());
    }

    if (h_.valid())
    {
        h_->reset();
    }
}


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

    // Get the mapper and the neighbouring patch
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& patchNbr = mapper.nbrFvPatch();

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

    if (qrName_ != "none")
    {
        sumq += patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    if (qrNbrName_ != "none")
    {
        sumq += mapper.fromNeighbour
        (
            patchNbr.lookupPatchField<volScalarField, scalar>(qrNbrName_)
        );
    }

    if (qrName_ != "none" || qrNbrName_ != "none")
    {
        sumq = qrRelax_*sumq + (1 - qrRelax_)*qrPrevious_;
        qrPrevious_ = sumq;
    }

    if (qs_.valid())
    {
        sumq += qs_();
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

        coupledTemperatureNbr.getNbr
        (
            sumKappaTByDeltaNbr,
            sumKappaByDeltaNbr,
            qNbr
        );

        // Include the effect of the optional neighbour insulation layer
        if (coupledTemperatureNbr.h_.valid())
        {
            const_cast<coupledTemperatureFvPatchScalarField&>
            (
                coupledTemperatureNbr
            ).h_->update();

            const scalarField hFactor
            (
                coupledTemperatureNbr.h_()
               /(coupledTemperatureNbr.h_() + sumKappaByDeltaNbr())
            );
            sumKappaTByDeltaNbr.ref() *= hFactor;
            sumKappaByDeltaNbr.ref() *= hFactor;
        }

        tmp<scalarField> sumKappaTByDeltaNbrMapped
        (
            mapper.fromNeighbour(sumKappaTByDeltaNbr)
        );

        tmp<scalarField> sumKappaByDeltaNbrMapped
        (
            mapper.fromNeighbour(sumKappaByDeltaNbr)
        );

        // Include the effect of the optional insulation layer
        if (h_.valid())
        {
            h_->update();

            const scalarField hFactor(h_()/(h_() + sumKappaByDeltaNbrMapped()));
            sumKappaTByDeltaNbrMapped.ref() *= hFactor;
            sumKappaByDeltaNbrMapped.ref() *= hFactor;
        }

        add(sumKappaTByDelta, sumKappaTByDeltaNbrMapped);
        add(sumKappaByDelta, sumKappaByDeltaNbrMapped);

        if (qNbr.valid())
        {
            sumq += mapper.fromNeighbour(qNbr);
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

        Info<< patch().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mapper.nbrMesh().name() << ':'
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

    if (qrName_ != "none")
    {
        writeEntry(os, "qr", qrName_);
        writeEntry(os, "qrRelaxation", qrRelax_);
        writeEntry(os, "qrPrevious", qrPrevious_);
    }

    if (Qs_ != 0)
    {
        writeEntry(os, "Qs", Qs_);
    }
    else if (qs_.valid())
    {
        writeEntry(os, "qs", qs_());
    }

    if (h_.valid())
    {
        writeEntry(os, h_());
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
