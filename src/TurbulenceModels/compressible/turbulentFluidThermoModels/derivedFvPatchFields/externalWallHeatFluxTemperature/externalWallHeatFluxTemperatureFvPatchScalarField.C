/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        externalWallHeatFluxTemperatureFvPatchScalarField::operationMode,
        3
    >::names[] =
    {
        "power",
        "flux",
        "coefficient"
    };
}

const Foam::NamedEnum
<
    Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationMode,
    3
> Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationModeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    mode_(fixedHeatFlux),
    Q_(0),
    Ta_(),
    relaxation_(1),
    emissivity_(0),
    qrRelaxation_(1),
    qrName_("undefined-qr"),
    thicknessLayers_(),
    kappaLayers_()
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.read(dict.lookup("mode"))),
    Q_(0),
    Ta_(),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1)),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0)),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    thicknessLayers_(),
    kappaLayers_()
{
    switch (mode_)
    {
        case fixedPower:
        {
            dict.lookup("Q") >> Q_;

            break;
        }
        case fixedHeatFlux:
        {
            q_ = scalarField("q", dict, p.size());

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_ = scalarField("h", dict, p.size());
            Ta_ = Function1<scalar>::New("Ta", dict);

            if (dict.found("thicknessLayers"))
            {
                dict.lookup("thicknessLayers") >> thicknessLayers_;
                dict.lookup("kappaLayers") >> kappaLayers_;
            }

            break;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (qrName_ != "none")
    {
        if (dict.found("qrPrevious"))
        {
            qrPrevious_ = scalarField("qrPrevious", dict, p.size());
        }
        else
        {
            qrPrevious_.setSize(p.size(), 0);
        }
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


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mode_(ptf.mode_),
    Q_(ptf.Q_),
    Ta_(ptf.Ta_, false),
    relaxation_(ptf.relaxation_),
    emissivity_(ptf.emissivity_),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_)
{
    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.setSize(mapper.size());
            q_.map(ptf.q_, mapper);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.setSize(mapper.size());
            h_.map(ptf.h_, mapper);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.setSize(mapper.size());
        qrPrevious_.map(ptf.qrPrevious_, mapper);
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    relaxation_(tppsf.relaxation_),
    emissivity_(tppsf.emissivity_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    relaxation_(tppsf.relaxation_),
    emissivity_(tppsf.emissivity_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.autoMap(m);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.autoMap(m);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.autoMap(m);
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const externalWallHeatFluxTemperatureFvPatchScalarField& tiptf =
        refCast<const externalWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.rmap(tiptf.q_, addr);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.rmap(tiptf.h_, addr);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.rmap(tiptf.qrPrevious_, addr);
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp(*this);

    scalarField qr(Tp.size(), 0);
    if (qrName_ != "none")
    {
        qr =
            qrRelaxation_
           *patch().lookupPatchField<volScalarField, scalar>(qrName_)
          + (1 - qrRelaxation_)*qrPrevious_;

        qrPrevious_ = qr;
    }

    switch (mode_)
    {
        case fixedPower:
        {
            refGrad() = (Q_/gSum(patch().magSf()) + qr)/kappa(Tp);
            refValue() = 0;
            valueFraction() = 0;

            break;
        }
        case fixedHeatFlux:
        {
            refGrad() = (q_ + qr)/kappa(Tp);
            refValue() = 0;
            valueFraction() = 0;

            break;
        }
        case fixedHeatTransferCoeff:
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
            scalarField hp(1/(1/h_ + totalSolidRes));

            const scalar Ta = Ta_->value(this->db().time().timeOutputValue());
            scalarField hpTa(hp*Ta);

            if (emissivity_ > 0)
            {
                // Evaluate the radiative flux to the environment
                // from the surface temperature ...
                if (totalSolidRes > 0)
                {
                    // ... including the effect of the solid wall thermal
                    // resistance
                    scalarField TpLambda(h_/(h_ + 1/totalSolidRes));
                    scalarField Ts(TpLambda*Tp + (1 - TpLambda)*Ta);
                    scalarField lambdaTa4(pow4((1 - TpLambda)*Ta));

                    hp += emissivity_*sigma.value()*(pow4(Ts) - lambdaTa4)/Tp;
                    hpTa += emissivity_*sigma.value()*(lambdaTa4 + pow4(Ta));
                }
                else
                {
                    // ... if there is no solid wall thermal resistance use
                    // the current wall temperature
                    hp += emissivity_*sigma.value()*pow3(Tp);
                    hpTa += emissivity_*sigma.value()*pow4(Ta);
                }
            }

            const scalarField kappaDeltaCoeffs
            (
                this->kappa(Tp)*patch().deltaCoeffs()
            );

            refGrad() = 0;

            forAll(Tp, i)
            {
                if (qr[i] < 0)
                {
                    const scalar hpmqr = hp[i] - qr[i]/Tp[i];

                    refValue()[i] = hpTa[i]/hpmqr;
                    valueFraction()[i] = hpmqr/(hpmqr + kappaDeltaCoeffs[i]);
                }
                else
                {
                    refValue()[i] = (hpTa[i] + qr[i])/hp[i];
                    valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
                }
            }

            break;
        }
    }

    valueFraction() = relaxation_*valueFraction() + (1 - relaxation_);
    refValue() = relaxation_*refValue() + (1 - relaxation_)*Tp;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

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


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    os.writeKeyword("mode")
        << operationModeNames[mode_] << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedPower:
        {
            os.writeKeyword("Q")
                << Q_ << token::END_STATEMENT << nl;

            break;
        }
        case fixedHeatFlux:
        {
            q_.writeEntry("q", os);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.writeEntry("h", os);
            Ta_->writeData(os);

            if (relaxation_ < 1)
            {
                os.writeKeyword("relaxation")
                    << relaxation_ << token::END_STATEMENT << nl;
            }

            if (emissivity_ > 0)
            {
                os.writeKeyword("emissivity")
                    << emissivity_ << token::END_STATEMENT << nl;
            }

            if (thicknessLayers_.size())
            {
                thicknessLayers_.writeEntry("thicknessLayers", os);
                kappaLayers_.writeEntry("kappaLayers", os);
            }

            break;
        }
    }

    os.writeKeyword("qr")<< qrName_ << token::END_STATEMENT << nl;

    if (qrName_ != "none")
    {
        os.writeKeyword("qrRelaxation")
            << qrRelaxation_ << token::END_STATEMENT << nl;

        qrPrevious_.writeEntry("qrPrevious", os);
    }

    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
