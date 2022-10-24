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

#include "externalTemperatureFvPatchScalarField.H"
#include "thermophysicalTransportModel.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalTemperatureFvPatchScalarField::
externalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    haveQ_(false),
    Q_(NaN),
    haveq_(false),
    q_(),
    haveh_(false),
    h_(),
    Ta_(),
    emissivity_(0),
    thicknessLayers_(),
    kappaLayers_(),
    relaxation_(1),
    qrName_(word::null),
    qrRelaxation_(1),
    qrPrevious_()
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


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
    Q_(haveQ_ ? dict.lookup<scalar>("Q") : NaN),
    haveq_(dict.found("q")),
    q_(haveq_ ? scalarField("q", dict, p.size()) : scalarField()),
    haveh_(dict.found("h")),
    h_(haveh_ ? scalarField("h", dict, p.size()) : scalarField()),
    Ta_(haveh_ ? Function1<scalar>::New("Ta", dict).ptr() : nullptr),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0)),
    thicknessLayers_
    (
        dict.lookupOrDefault<scalarList>("thicknessLayers", scalarList())
    ),
    kappaLayers_
    (
        dict.lookupOrDefault<scalarList>("kappaLayers", scalarList())
    ),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", word::null)),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
    qrPrevious_
    (
        qrName_ != word::null
      ? dict.found("qrPrevious")
      ? scalarField("qrPrevious", dict, p.size())
      : scalarField(0, p.size())
      : scalarField()
    )
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (!haveQ_ && !haveq_ && !haveh_)
    {
        FatalIOErrorInFunction(dict)
            << "One or more of Q (heat power), q (heat flux), and h (heat "
            << "transfer coefficient) must be specified"
            << exit(FatalIOError);
    }

    if (thicknessLayers_.size() != kappaLayers_.size())
    {
        FatalIOErrorInFunction(dict)
            << "If either thicknessLayers or kappaLayers is specified, then "
            << "both must be specified and be lists of the same length "
            << exit(FatalIOError);
    }

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


Foam::externalTemperatureFvPatchScalarField::
externalTemperatureFvPatchScalarField
(
    const externalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    haveQ_(ptf.haveQ_),
    Q_(ptf.Q_),
    haveq_(ptf.haveq_),
    q_(haveq_ ? mapper(ptf.q_)() : scalarField()),
    haveh_(ptf.haveh_),
    h_(haveh_ ? mapper(ptf.h_)() : scalarField()),
    Ta_(ptf.Ta_, false),
    emissivity_(ptf.emissivity_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    relaxation_(ptf.relaxation_),
    qrName_(ptf.qrName_),
    qrRelaxation_(ptf.qrRelaxation_),
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
    Q_(tppsf.Q_),
    haveq_(tppsf.haveq_),
    q_(tppsf.q_),
    haveh_(tppsf.haveh_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    emissivity_(tppsf.emissivity_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
    relaxation_(tppsf.relaxation_),
    qrName_(tppsf.qrName_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrPrevious_(tppsf.qrPrevious_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    if (haveq_)
    {
        m(q_, q_);
    }

    if (haveh_)
    {
        m(h_, h_);
    }

    if (qrName_ != word::null)
    {
        m(qrPrevious_, qrPrevious_);
    }
}


void Foam::externalTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const externalTemperatureFvPatchScalarField& tiptf =
        refCast<const externalTemperatureFvPatchScalarField>(ptf);

    if (haveq_)
    {
        q_.rmap(tiptf.q_, addr);
    }

    if (haveh_)
    {
        h_.rmap(tiptf.h_, addr);
    }

    if (qrName_ != word::null)
    {
        qrPrevious_.rmap(tiptf.qrPrevious_, addr);
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

    if (haveq_)
    {
        q_.reset(tiptf.q_);
    }

    if (haveh_)
    {
        h_.reset(tiptf.h_);
    }

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

    const scalarField& Tp(*this);

    // Store current valueFraction and refValue for relaxation
    const scalarField valueFraction0(valueFraction());
    const scalarField refValue0(refValue());

    // Get the radiative heat flux and relax
    scalarField qr(Tp.size(), 0);
    if (qrName_ != word::null)
    {
        qr =
            qrRelaxation_
           *patch().lookupPatchField<volScalarField, scalar>(qrName_)
          + (1 - qrRelaxation_)*qrPrevious_;
        qrPrevious_ = qr;
    }

    // Compute the total non-convective heat flux
    scalarField qTot(qr);
    if (haveQ_)
    {
        qTot += Q_/gSum(patch().magSf());
    }
    if (haveq_)
    {
        qTot += q_;
    }

    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    const scalarField kappa(ttm.kappaEff(patch().index()));
    tmp<scalarField> qCorr(ttm.qCorr(patch().index()));

    if (qCorr.valid())
    {
        qTot += qCorr;
    }

    // Evaluate
    if (!haveh_)
    {
        refGrad() = qTot/kappa;
        refValue() = Tp;
        valueFraction() = 0;
    }
    else
    {
        scalar totalSolidRes = 0;
        if (thicknessLayers_.size())
        {
            forAll(thicknessLayers_, iLayer)
            {
                const scalar l = thicknessLayers_[iLayer];
                if (kappaLayers_[iLayer] > 0)
                {
                    totalSolidRes += l/kappaLayers_[iLayer];
                }
            }
        }

        const scalar Ta = Ta_->value(this->db().time().userTimeValue());

        const scalarField hp
        (
            1
           /(
                1
               /(
                    (emissivity_ > 0)
                  ? (
                        h_
                      + emissivity_*sigma.value()
                       *((pow3(Ta) + pow3(Tp)) + Ta*Tp*(Ta + Tp))
                    )()
                  : h_
                ) + totalSolidRes
            )
        );

        const scalarField hpTa(hp*Ta);

        const scalarField kappaDeltaCoeffs
        (
            kappa*patch().deltaCoeffs()
        );

        refGrad() = 0;
        forAll(Tp, i)
        {
            if (qTot[i] < 0)
            {
                const scalar hpmqTot = hp[i] - qTot[i]/Tp[i];
                refValue()[i] = hpTa[i]/hpmqTot;
                valueFraction()[i] = hpmqTot/(hpmqTot + kappaDeltaCoeffs[i]);
            }
            else
            {
                refValue()[i] = (hpTa[i] + qTot[i])/hp[i];
                valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
            }
        }
    }

    // Relax
    valueFraction() =
        relaxation_*valueFraction() + (1 - relaxation_)*valueFraction0;
    refValue() =
        relaxation_*refValue() + (1 - relaxation_)*refValue0;

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
        writeEntry(os, "Q", Q_);
    }

    if (haveq_)
    {
        writeEntry(os, "q", q_);
    }

    if (haveh_)
    {
        writeEntry(os, "h", h_);
        writeEntry(os, Ta_());
        writeEntryIfDifferent(os, "emissivity", scalar(0), emissivity_);
        writeEntryIfDifferent
        (
            os,
            "thicknessLayers",
            scalarList(),
            thicknessLayers_
        );
        writeEntryIfDifferent
        (
            os,
            "kappaLayers",
            scalarList(),
            kappaLayers_
        );
    }

    writeEntryIfDifferent(os, "relaxation", scalar(1), relaxation_);

    if (qrName_ != word::null)
    {
        writeEntry(os, "qr", qrName_);
        writeEntry(os, "qrRelaxation", qrRelaxation_);
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
        patch,
        externalWallHeatFluxTemperature,
        "externalWallHeatFluxTemperature"
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
