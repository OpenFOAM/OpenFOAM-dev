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

#include "thermalBaffle1DFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "thermophysicalTransportModel.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    qs_(p.size()),
    solidDict_(),
    solidPtr_(nullptr),
    qrPrevious_(p.size()),
    qrRelaxation_(1),
    qrName_("undefined-qr")
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    TName_("T"),
    baffleActivated_(dict.lookupOrDefault<bool>("baffleActivated", true)),
    thickness_(),
    qs_(p.size(), 0),
    solidDict_(dict),
    solidPtr_(),
    qrPrevious_(p.size(), 0.0),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", "none"))
{
    mappedPatchBase::validateMapForField
    (
        *this,
        dict,
        mappedPatchBase::from::sameRegion
      & mappedPatchBase::from::differentPatch
    );

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dict, p.size());
    }

    if (dict.found("qs"))
    {
        qs_ = scalarField("qs", dict, p.size());
    }
    else if (dict.found("Qs"))
    {
        qs_ = scalarField("Qs", dict, p.size());
    }

    if (dict.found("qrPrevious"))
    {
        qrPrevious_ = scalarField("qrPrevious", dict, p.size());
    }

    if (dict.found("refValue") && baffleActivated_)
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }
}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(mapper(ptf.thickness_)),
    qs_(mapper(ptf.qs_)),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    qrPrevious_(mapper(ptf.qrPrevious_)),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    qrPrevious_(ptf.qrPrevious_),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class solidType>
bool thermalBaffle1DFvPatchScalarField<solidType>::owner() const
{
    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    return patch().patch().index() < mpp.nbrPolyPatch().index();
}


template<class solidType>
const thermalBaffle1DFvPatchScalarField<solidType>&
thermalBaffle1DFvPatchScalarField<solidType>::nbrField() const
{
    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const polyMesh& nbrMesh = mpp.nbrMesh();
    const label nbrPatchi = mpp.nbrPolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchi];

    return
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>
            (
                TName_
            )
        );
}


template<class solidType>
const solidType& thermalBaffle1DFvPatchScalarField<solidType>::solid() const
{
    if (this->owner())
    {
        if (solidPtr_.empty())
        {
            solidPtr_.reset(new solidType(solidDict_));
        }

        return solidPtr_();
    }
    else
    {
        return nbrField().solid();
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::
baffleThickness() const
{
    if (this->owner())
    {
        if (thickness_.size() != patch().size())
        {
            FatalIOErrorInFunction(solidDict_)
                << " Field thickness has not been specified "
                << " for patch " << this->patch().name()
                << exit(FatalIOError);
        }

        return thickness_;
    }
    else
    {
        const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
        return mpp.distribute(nbrField().baffleThickness());
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::qs() const
{
    if (this->owner())
    {
         return qs_;
    }
    else
    {
        const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
        return mpp.distribute(nbrField().qs());
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    if (this->owner())
    {
        m(thickness_, thickness_);
        m(qs_, qs_);
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.rmap(tiptf.thickness_, addr);
        qs_.rmap(tiptf.qs_, addr);
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.reset(tiptf.thickness_);
        qs_.reset(tiptf.qs_);
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());

    if (baffleActivated_)
    {
        const thermophysicalTransportModel& ttm =
            db().objectRegistry::template lookupObject
            <
                thermophysicalTransportModel
            >
            (
                IOobject::groupName
                (
                    thermophysicalTransportModel::typeName,
                    internalField().group()
                )
            );

        // Local properties
        const fvPatchScalarField& Tp =
            patch().template lookupPatchField<volScalarField, scalar>(TName_);

        const scalarField kappap(ttm.kappaEff(patch().index()));

        scalarField qr(Tp.size(), Zero);

        if (qrName_ != "none")
        {
            qr = patch().template lookupPatchField<volScalarField, scalar>
                (qrName_);

            qr = qrRelaxation_*qr + (1.0 - qrRelaxation_)*qrPrevious_;
            qrPrevious_ = qr;
        }

        const scalarField kappaDelta(kappap*patch().deltaCoeffs());

        // Neighbour properties
        const scalarField nbrTp(mpp.distribute(nbrField()));

        // Solid properties
        scalarField kappas(patch().size(), 0.0);
        forAll(kappas, i)
        {
            kappas[i] = solid().kappa(0.0, (Tp[i] + nbrTp[i])/2.0);
        }

        const scalarField kappaByDeltaSolid(kappas/baffleThickness());

        const scalarField alpha(kappaByDeltaSolid - qr/Tp);

        valueFraction() = alpha/(alpha + kappaDelta);

        refValue() = (kappaByDeltaSolid*nbrTp + qs()/2.0)/alpha;

        if (debug)
        {
            scalar Q = gAverage(kappap*snGrad());
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrField().patch().name() << ':'
                << this->internalField().name() << " :"
                << " heat[W]:" << Q
                << " walltemperature "
                << " min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << endl;
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);

    if (this->owner())
    {
        writeEntry(os, "thickness", baffleThickness()());
        writeEntry(os, "qs", qs()());
        solid().write(os);
    }

    writeEntry(os, "qrPrevious", qrPrevious_);
    writeEntry(os, "qr", qrName_);
    writeEntry(os, "qrRelaxation", qrRelaxation_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
