/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"
#include "thermophysicalTransportModel.H"
#include "mapDistribute.H"

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
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    Qs_(p.size()),
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
    const thermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase(p.patch(), ptf),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(mapper(ptf.thickness_)),
    Qs_(mapper(ptf.Qs_)),
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
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), NEARESTPATCHFACE, dict),
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(dict.lookupOrDefault<bool>("baffleActivated", true)),
    thickness_(),
    Qs_(p.size(), 0),
    solidDict_(dict),
    solidPtr_(),
    qrPrevious_(p.size(), 0.0),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", "none"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dict, p.size());
    }

    if (dict.found("Qs"))
    {
        Qs_ = scalarField("Qs", dict, p.size());
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
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    Qs_(ptf.Qs_),
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
    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    return (patchi < nbrPatchi);
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
        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        return nbrField.solid();
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
            FatalIOErrorInFunction
            (
                solidDict_
            )<< " Field thickness has not been specified "
            << " for patch " << this->patch().name()
            << exit(FatalIOError);
        }

        return thickness_;
    }
    else
    {
        const mapDistribute& mapDist = this->mappedPatchBase::map();

        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];
        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        tmp<scalarField> tthickness
        (
            new scalarField(nbrField.baffleThickness())
        );
        scalarField& thickness = tthickness.ref();
        mapDist.distribute(thickness);
        return tthickness;
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::Qs() const
{
    if (this->owner())
    {
         return Qs_;
    }
    else
    {
        const mapDistribute& mapDist = this->mappedPatchBase::map();

        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        tmp<scalarField> tQs(new scalarField(nbrField.Qs()));
        scalarField& Qs = tQs.ref();
        mapDist.distribute(Qs);
        return tQs;
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mappedPatchBase::clearOut();

    mixedFvPatchScalarField::autoMap(m);

    if (this->owner())
    {
        m(thickness_, thickness_);
        m(Qs_, Qs_);
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
        Qs_.rmap(tiptf.Qs_, addr);
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
    UPstream::msgType() = oldTag+1;

    const mapDistribute& mapDist = this->mappedPatchBase::map();

    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    if (baffleActivated_)
    {
        const fvPatch& nbrPatch = patch().boundaryMesh()[nbrPatchi];

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

        // local properties
        const scalarField kappaw(ttm.kappaEff(patchi));

        const fvPatchScalarField& Tp =
            patch().template lookupPatchField<volScalarField, scalar>(TName_);


        scalarField qr(Tp.size(), 0.0);

        if (qrName_ != "none")
        {
            qr = patch().template lookupPatchField<volScalarField, scalar>
                (qrName_);

            qr = qrRelaxation_*qr + (1.0 - qrRelaxation_)*qrPrevious_;
            qrPrevious_ = qr;
        }

        tmp<scalarField> Ti = patchInternalField();

        scalarField myKDelta(patch().deltaCoeffs()*kappaw);

        // nrb properties
        scalarField nbrTp = ttm.thermo().T().boundaryField()[nbrPatchi];
        mapDist.distribute(nbrTp);

        // solid properties
        scalarField kappas(patch().size(), 0.0);
        forAll(kappas, i)
        {
            kappas[i] = solid().kappa(0.0, (Tp[i] + nbrTp[i])/2.0);
        }

        scalarField KDeltaSolid(kappas/baffleThickness());

        scalarField alpha(KDeltaSolid - qr/Tp);

        valueFraction() = alpha/(alpha + myKDelta);

        refValue() = (KDeltaSolid*nbrTp + Qs()/2.0)/alpha;

        if (debug)
        {
            scalar Q = gAverage(kappaw*snGrad());
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrPatch.name() << ':'
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
    mappedPatchBase::write(os);

    if (this->owner())
    {
        writeEntry(os, "thickness", baffleThickness()());
        writeEntry(os, "Qs", Qs()());
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
