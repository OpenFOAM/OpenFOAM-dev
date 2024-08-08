/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "mappedFvPatchBaseBase.H"

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
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    TName_("T"),
    baffleActivated_(dict.lookupOrDefault<bool>("baffleActivated", true)),
    thickness_(),
    qs_(p.size(), 0),
    solidPtr_
    (
        this->owner()
      ? new solidType("solid", dict)
      : nullptr
    ),
    qrPrevious_(p.size(), 0.0),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", dimless, 1)),
    qrName_(dict.lookupOrDefault<word>("qr", "none"))
{
    mappedPatchBaseBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBaseBase::from::sameRegion
      & mappedPatchBaseBase::from::differentPatch
    );

    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dimLength, dict, p.size());
    }

    if (dict.found("qs"))
    {
        qs_ = scalarField("qs", dimPower/dimArea, dict, p.size());
    }

    if (dict.found("qrPrevious"))
    {
        qrPrevious_ =
            scalarField("qrPrevious", dimPower/dimArea, dict, p.size());
    }

    if (dict.found("refValue") && baffleActivated_)
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
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(mapper(ptf.thickness_)),
    qs_(mapper(ptf.qs_)),
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
    solidPtr_(ptf.solidPtr_),
    qrPrevious_(ptf.qrPrevious_),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class solidType>
bool thermalBaffle1DFvPatchScalarField<solidType>::owner() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    return patch().index() < mapper.nbrFvPatch().index();
}


template<class solidType>
const thermalBaffle1DFvPatchScalarField<solidType>&
thermalBaffle1DFvPatchScalarField<solidType>::nbrField() const
{
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const polyMesh& nbrMesh = mapper.nbrMesh();
    const label nbrPatchi = mapper.nbrFvPatch().index();
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
            FatalErrorInFunction
                << "solid not allocated" << exit(FatalError);
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
            FatalErrorInFunction
                << " Field thickness has not been specified "
                << " for patch " << this->patch().name()
                << exit(FatalError);
        }

        return thickness_;
    }
    else
    {
        const mappedFvPatchBaseBase& mapper =
            mappedFvPatchBaseBase::getMap(patch());
        return mapper.fromNeighbour(nbrField().baffleThickness());
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
        const mappedFvPatchBaseBase& mapper =
            mappedFvPatchBaseBase::getMap(patch());
        return mapper.fromNeighbour(nbrField().qs());
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        mapper(thickness_, tiptf.thickness_);
        mapper(qs_, tiptf.qs_);
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

    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());

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
        const scalarField nbrTp(mapper.fromNeighbour(nbrField()));

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
