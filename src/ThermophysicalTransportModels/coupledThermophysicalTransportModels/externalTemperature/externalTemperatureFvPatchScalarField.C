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

#include "externalTemperatureFvPatchScalarField.H"
#include "thermophysicalTransportModel.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::externalTemperatureFvPatchScalarField::getKappa
(
    scalarField& kappa,
    scalarField& sumKappaTByDelta,
    scalarField& sumKappaByDelta,
    scalarField& Tref,
    scalarField& Tw,
    scalarField& sumq,
    scalarField& qByKappa
) const
{
    const thermophysicalTransportModel& ttm =
        patch().boundaryMesh().mesh()
       .lookupType<thermophysicalTransportModel>();

    kappa = ttm.kappaEff(patch().index());

    tmp<scalarField> qCorr(ttm.qCorr(patch().index()));

    if (qCorr.valid())
    {
        sumq += qCorr;
    }

    qByKappa = sumq/kappa;

    Tw = *this;
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
      ? Function1<scalar>::New("Q", dict)
      : autoPtr<Function1<scalar>>()
    ),
    haveq_(dict.found("q")),
    q_
    (
        haveq_
      ? Function1<scalar>::New("q", dict)
      : autoPtr<Function1<scalar>>()
    ),
    haveh_(dict.found("h")),
    h_(haveh_ ? Function1<scalar>::New("h", dict).ptr() : nullptr),
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
    relax_(dict.lookupOrDefault<scalar>("relaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", word::null)),
    qrRelax_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
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
    Q_(ptf.Q_, false),
    haveq_(ptf.haveq_),
    q_(ptf.q_, false),
    haveh_(ptf.haveh_),
    h_(ptf.h_, false),
    Ta_(ptf.Ta_, false),
    emissivity_(ptf.emissivity_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
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
    Ta_(tppsf.Ta_, false),
    emissivity_(tppsf.emissivity_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
    relax_(tppsf.relax_),
    qrName_(tppsf.qrName_),
    qrRelax_(tppsf.qrRelax_),
    qrPrevious_(tppsf.qrPrevious_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
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

    // Get the radiative heat flux and relax
    scalarField qr(size(), 0);
    if (qrName_ != word::null)
    {
        qr =
            qrRelax_
           *patch().lookupPatchField<volScalarField, scalar>(qrName_)
          + (1 - qrRelax_)*qrPrevious_;
        qrPrevious_ = qr;
    }

    // Compute the total non-convective heat flux
    scalarField sumq(qr);
    if (haveQ_)
    {
        sumq += Q_->value(db().time().userTimeValue())/gSum(patch().magSf());
    }
    if (haveq_)
    {
        sumq += q_->value(db().time().userTimeValue());
    }

    scalarField kappa(size(), 0);
    scalarField sumKappaTByDelta(size(), 0);
    scalarField sumKappaByDelta(size(), 0);
    scalarField Tref(*this);
    scalarField Tw(*this);
    scalarField qByKappa(size(), 0);
    getKappa
    (
        kappa,
        sumKappaTByDelta,
        sumKappaByDelta,
        Tref,
        Tw,
        sumq,
        qByKappa
    );

    // Add optional external convective heat transfer contribution
    if (haveh_)
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

        const scalar h = h_->value(this->db().time().userTimeValue());
        const scalar Ta = Ta_->value(this->db().time().userTimeValue());

        const scalarField hp
        (
            1
           /(
                1
               /(
                    (emissivity_ > 0)
                  ? (
                        h
                      + emissivity_*sigma.value()
                       *((pow3(Ta) + pow3(Tw)) + Ta*Tw*(Ta + Tw))
                    )()
                  : scalarField(size(), h)
                ) + totalSolidRes
            )
        );

        sumKappaByDelta += hp;
        sumKappaTByDelta += hp*Ta;

        refValue() = sumKappaTByDelta/sumKappaByDelta;
    }
    else
    {
        refValue() = Tref;
    }

    valueFraction() =
        sumKappaByDelta/(kappa*patch().deltaCoeffs() + sumKappaByDelta);

    refGrad() = qByKappa;

    // Modify mixed parameters for under-relaxation
    if (relax_ != 1)
    {
        const scalarField f(valueFraction());
        valueFraction() = 1 - relax_*(1 - f);
        refValue() = (f*relax_*refValue() + (1 - relax_)*Tw)/valueFraction();
        // refGrad() = No change
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
        writeEntry(os, Q_());
    }

    if (haveq_)
    {
        writeEntry(os, q_());
    }

    if (haveh_)
    {
        writeEntry(os, h_());
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
