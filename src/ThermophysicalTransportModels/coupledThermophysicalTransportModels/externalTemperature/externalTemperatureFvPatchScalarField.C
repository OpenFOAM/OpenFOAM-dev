/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "externalTemperatureFvPatchScalarField.H"
#include "thermophysicalTransportModel.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::externalTemperatureFvPatchScalarField::plusEqOp
(
    tmp<scalarField>& tf,
    const scalar d
) const
{
    if (!tf.valid())
    {
        tf = new scalarField(size(), d);
    }
    else
    {
        tf.ref() += d;
    }
}


void Foam::externalTemperatureFvPatchScalarField::plusEqOp
(
    tmp<scalarField>& tf,
    const tmp<scalarField>& tdf
) const
{
    if (!tdf.valid())
    {
        return;
    }

    if (!tf.valid())
    {
        tf = tdf.ptr();
    }
    else
    {
        tf.ref() += tdf;
    }
}


void Foam::externalTemperatureFvPatchScalarField::getKappa
(
    scalarField& kappa,
    tmp<scalarField>& sumKappaTcByDelta,
    tmp<scalarField>& sumKappaByDelta,
    tmp<scalarField>& T,
    tmp<scalarField>& sumq
) const
{
    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    kappa = ttm.kappaEff(patch().index());

    T = tmp<scalarField>(*this);

    plusEqOp(sumq, ttm.qCorr(patch().index()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalTemperatureFvPatchScalarField::
externalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    haveQ_(dict.found("Q")),
    Q_
    (
        haveQ_
      ? Function1<scalar>::New("Q", db().time().userUnits(), dimPower, dict)
      : autoPtr<Function1<scalar>>()
    ),
    haveq_(dict.found("q")),
    q_
    (
        haveq_
      ? Function1<scalar>::New
        (
            "q",
            db().time().userUnits(),
            dimPower/dimArea,
            dict
        )
      : autoPtr<Function1<scalar>>()
    ),
    haveh_(dict.found("h")),
    h_
    (
        haveh_
      ? Function1<scalar>::New
        (
            "h",
            db().time().userUnits(),
            dimPower/dimArea/dimTemperature,
            dict
        ).ptr()
      : nullptr
    ),
    haveEmissivity_(dict.found("emissivity")),
    emissivity_
    (
        haveEmissivity_
      ? dict.lookup<scalar>("emissivity", unitFraction)
      : NaN
    ),
    haveLayers_(dict.found("thicknessLayers") || dict.found("kappaLayers")),
    thicknessLayers_
    (
        haveLayers_
      ? dict.lookup<scalarList>("thicknessLayers", dimLength)
      : scalarList()
    ),
    kappaLayers_
    (
        haveLayers_
      ? dict.lookup<scalarList>("kappaLayers", dimThermalConductivity)
      : scalarList()
    ),
    Ta_
    (
        haveh_ || haveEmissivity_
      ? Function1<scalar>::New
        (
            "Ta",
            db().time().userUnits(),
            dimTemperature,
            dict
        ).ptr()
      : nullptr
    ),
    relax_(dict.lookupOrDefault<scalar>("relaxation", unitFraction, 1)),
    qrName_(dict.lookupOrDefault<word>("qr", word::null)),
    qrRelax_(dict.lookupOrDefault<scalar>("qrRelaxation", unitFraction, 1)),
    qrPrevious_
    (
        qrName_ != word::null
      ? dict.found("qrPrevious")
      ? scalarField("qrPrevious", dimPower/dimArea, dict, p.size())
      : scalarField(p.size(), 0)
      : scalarField()
    )
{
    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    if (haveEmissivity_ && (emissivity_ < 0 || emissivity_ > 1))
    {
        FatalIOErrorInFunction(dict)
            << "Emissivity must be in the range 0 to 1"
            << exit(FatalIOError);
    }

    if (thicknessLayers_.size() != kappaLayers_.size())
    {
        FatalIOErrorInFunction(dict)
            << "If either thicknessLayers or kappaLayers is specified, then "
            << "both must be specified and be lists of the same length "
            << exit(FatalIOError);
    }

    if (haveEmissivity_ && haveLayers_)
    {
        FatalIOErrorInFunction(dict)
            << "Emissivity and thicknessLayers/kappaLayers are incompatible"
            << exit(FatalIOError);
    }

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", iF.dimensions(), dict, p.size());
        refGrad() =
            scalarField
            (
                "refGradient",
                iF.dimensions()/dimLength,
                dict,
                p.size()
            );
        valueFraction() =
            scalarField("valueFraction", unitFraction, dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::externalTemperatureFvPatchScalarField::
externalTemperatureFvPatchScalarField
(
    const externalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    haveQ_(ptf.haveQ_),
    Q_(ptf.Q_, false),
    haveq_(ptf.haveq_),
    q_(ptf.q_, false),
    haveh_(ptf.haveh_),
    h_(ptf.h_, false),
    haveEmissivity_(ptf.haveEmissivity_),
    emissivity_(ptf.emissivity_),
    haveLayers_(ptf.haveLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    Ta_(ptf.Ta_, false),
    relax_(ptf.relax_),
    qrName_(ptf.qrName_),
    qrRelax_(ptf.qrRelax_),
    qrPrevious_
    (
        qrName_ != word::null
      ? mapper(ptf.qrPrevious_)()
      : scalarField()
    )
{}


Foam::externalTemperatureFvPatchScalarField::
externalTemperatureFvPatchScalarField
(
    const externalTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    haveQ_(tppsf.haveQ_),
    Q_(tppsf.Q_, false),
    haveq_(tppsf.haveq_),
    q_(tppsf.q_, false),
    haveh_(tppsf.haveh_),
    h_(tppsf.h_, false),
    haveEmissivity_(tppsf.haveEmissivity_),
    emissivity_(tppsf.emissivity_),
    haveLayers_(tppsf.haveLayers_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
    Ta_(tppsf.Ta_, false),
    relax_(tppsf.relax_),
    qrName_(tppsf.qrName_),
    qrRelax_(tppsf.qrRelax_),
    qrPrevious_(tppsf.qrPrevious_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);

    const externalTemperatureFvPatchScalarField& tiptf =
        refCast<const externalTemperatureFvPatchScalarField>(ptf);

    if (qrName_ != word::null)
    {
        mapper(qrPrevious_, tiptf.qrPrevious_);
    }
}


void Foam::externalTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const externalTemperatureFvPatchScalarField& tiptf =
        refCast<const externalTemperatureFvPatchScalarField>(ptf);

    if (qrName_ != word::null)
    {
        qrPrevious_.reset(tiptf.qrPrevious_);
    }
}



void Foam::externalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().value();

    // Get thermal conductivities, the patch temperature and any explicit
    // correction heat flux
    scalarField kappa;
    tmp<scalarField> sumKappaTcByDelta, sumKappaByDelta;
    tmp<scalarField> T;
    tmp<scalarField> sumq;
    getKappa(kappa, sumKappaTcByDelta, sumKappaByDelta, T, sumq);

    // Add any user specified heat fluxes
    if (haveQ_)
    {
        plusEqOp(sumq, Q_->value(t)/gSum(patch().magSf()));
    }
    if (haveq_)
    {
        plusEqOp(sumq, q_->value(t));
    }

    // Add the (relaxed) radiative heat flux
    if (qrName_ != word::null)
    {
        const fvPatchScalarField& qrCurrent =
            patch().lookupPatchField<volScalarField, scalar>(qrName_);

        const scalarField qr(qrRelax_*qrCurrent + (1 - qrRelax_)*qrPrevious_);

        qrPrevious_ = qr;

        plusEqOp(sumq, qr);
    }

    // Evaluate the ambient temperature
    const scalar Ta = haveh_ || haveEmissivity_ ? Ta_->value(t) : NaN;

    // Evaluate the combined convective and radiative heat transfer coefficient
    tmp<scalarField> hEff;
    if (haveh_)
    {
        plusEqOp(hEff, h_->value(t));
    }
    if (haveEmissivity_)
    {
        plusEqOp
        (
            hEff,
            emissivity_
           *constant::physicoChemical::sigma.value()
           *(sqr(Ta) + sqr(T()))
           *(Ta + T())
        );
    }

    // Determine the (reciprocal of the) heat transfer coefficient for the
    // layer resistances and combine with the convective and radiative heat
    // transfer coefficients to create a complete effective coefficient
    if (hEff.valid() && haveLayers_)
    {
        scalar oneByHLayers = 0;
        if (haveLayers_)
        {
            forAll(thicknessLayers_, layeri)
            {
                oneByHLayers += thicknessLayers_[layeri]/kappaLayers_[layeri];
            }
        }

        hEff = 1/(1/hEff + oneByHLayers);
    }

    // If we have a heat transfer coefficient then add it to the kappa sums
    if (hEff.valid())
    {
        plusEqOp(sumKappaByDelta, hEff());
        plusEqOp(sumKappaTcByDelta, hEff*Ta);
    }

    // Set the mixed parameters
    const scalarField kappaByDelta(kappa*patch().deltaCoeffs());
    tmp<scalarField> kappaPlusSumKappaByDelta
    (
        sumKappaByDelta.valid()
      ? kappaByDelta + sumKappaByDelta()
      : tmp<scalarField>(kappaByDelta)
    );

    // ... value fraction
    if (sumKappaByDelta.valid())
    {
        valueFraction() = sumKappaByDelta()/kappaPlusSumKappaByDelta();
    }
    else
    {
        valueFraction() = Zero;
    }

    // ... reference value
    tmp<scalarField> trefValue;
    if (sumKappaByDelta.valid())
    {
        plusEqOp
        (
            trefValue,
            max(sumKappaTcByDelta, small*kappaByDelta*internalField())
           /max(sumKappaByDelta, small*kappaByDelta)
        );
    }
    if (sumq.valid())
    {
        plusEqOp(trefValue, sumq()/kappaPlusSumKappaByDelta());
    }
    if (trefValue.valid())
    {
        refValue() = trefValue;
    }
    else
    {
        refValue() = Zero;
    }

    // ... and reference gradient
    if (sumq.valid())
    {
        refGrad() = sumq*patch().deltaCoeffs()/kappaPlusSumKappaByDelta();
    }
    else
    {
        refGrad() = Zero;
    }

    // Modify the mixed parameters for under-relaxation
    if (relax_ != 1)
    {
        const scalarField f(valueFraction());
        valueFraction() = 1 - relax_*(1 - f);
        refValue() = (f*relax_*refValue() + (1 - relax_)*T)/valueFraction();
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::externalTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    if (haveQ_)
    {
        writeEntry(os, db().time().userUnits(), dimPower, Q_());
    }

    if (haveq_)
    {
        writeEntry(os, db().time().userUnits(), dimPower/dimArea, q_());
    }

    if (haveh_)
    {
        writeEntry
        (
            os,
            db().time().userUnits(),
            dimPower/dimArea/dimTemperature,
            h_()
        );
    }

    if (haveEmissivity_)
    {
        writeEntry(os, "emissivity", emissivity_);
    }

    if (haveLayers_)
    {
        writeEntry(os, "thicknessLayers", thicknessLayers_);
        writeEntry(os, "kappaLayers", kappaLayers_);
    }

    if (haveh_ || haveEmissivity_)
    {
        writeEntry(os, db().time().userUnits(), dimTemperature, Ta_());
    }

    writeEntryIfDifferent(os, "relaxation", scalar(1), relax_);

    if (qrName_ != word::null)
    {
        writeEntry(os, "qr", qrName_);
        writeEntry(os, "qrRelaxation", qrRelax_);
        writeEntry(os, "qrPrevious", qrPrevious_);
    }

    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalTemperatureFvPatchScalarField
    );

    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        externalTemperatureFvPatchScalarField,
        patchMapper,
        externalWallHeatFluxTemperature,
        "externalWallHeatFluxTemperature"
    );

    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvPatchScalarField,
        externalTemperatureFvPatchScalarField,
        dictionary,
        externalWallHeatFluxTemperature,
        "externalWallHeatFluxTemperature"
    );
}


// ************************************************************************* //
