/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "energyRegionCoupledFvPatchScalarField.H"
#include "Time.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::energyRegionCoupledFvPatchScalarField::kappaMethodType,
        3
    >::names[] =
    {
        "solid",
        "fluid",
        "undefined"
    };
}


const Foam::NamedEnum
<
    Foam::energyRegionCoupledFvPatchScalarField::kappaMethodType,
    3
> Foam::energyRegionCoupledFvPatchScalarField::methodTypeNames_;


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//

void Foam::energyRegionCoupledFvPatchScalarField::setMethod() const
{
    if (method_ == UNDEFINED)
    {
        if
        (
            this->db().foundObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            )
        )
        {
            method_ = FLUID;
        }
        else
        {
            method_ = SOLID;
        }
    }

    if (!nbrThermoPtr_)
    {
        nbrThermoPtr_ =
        (
            &regionCoupledPatch_.nbrMesh().lookupObject<basicThermo>
            (
                basicThermo::dictName
            )
        );
    }

    if (!thermoPtr_)
    {
        thermoPtr_ =
        (
            &this->db().lookupObject<basicThermo>
            (
                basicThermo::dictName
            )
        );
    }
}


Foam::tmp<Foam::scalarField> Foam::energyRegionCoupledFvPatchScalarField::
kappa() const
{
    switch (method_)
    {
        case FLUID:
        {
            const compressible::turbulenceModel& turbModel =
                this->db().lookupObject<compressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            return turbModel.kappaEff(patch().index());
        }
        break;

        case SOLID:
        {
            const basicThermo& thermo =
                this->db().lookupObject<basicThermo>
                (
                    basicThermo::dictName
                );

            return thermo.kappa(patch().index());
        }
        break;

        case UNDEFINED:
        {
            FatalErrorInFunction
                    << " on mesh " << this->db().name() << " patch "
                    << patch().name()
                    << " could not find a method in. Methods are:  "
                    << methodTypeNames_.toc()
                    << " Not turbulenceModel or thermophysicalProperties"
                    << " were found"
                    << exit(FatalError);
        }
        break;
    }
    return scalarField(0);
}


Foam::tmp<Foam::scalarField> Foam::energyRegionCoupledFvPatchScalarField::
weights() const
{
    const fvPatch& patch = regionCoupledPatch_.patch();

    const scalarField deltas
    (
        patch.nf() & patch.delta()
    );

    const scalarField alphaDelta(kappa()/deltas);

    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

    const energyRegionCoupledFvPatchScalarField& nbrField =
    refCast
    <
        const energyRegionCoupledFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>("T")
    );

    // Needed in the first calculation of weights
    nbrField.setMethod();

    const scalarField nbrAlpha
    (
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrField.kappa()
        )
    );

    const scalarField nbrDeltas
    (
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrPatch.nf() & nbrPatch.delta()
        )
    );

    const scalarField nbrAlphaDelta(nbrAlpha/nbrDeltas);

    tmp<scalarField> tw(new scalarField(deltas.size()));
    scalarField& w = tw.ref();

    forAll(alphaDelta, facei)
    {
        scalar di = alphaDelta[facei];
        scalar dni = nbrAlphaDelta[facei];

        w[facei] = di/(di + dni);
    }

    return tw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyRegionCoupledFvPatchScalarField::
energyRegionCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledFvPatchField<scalar>(p, iF),
    regionCoupledPatch_(refCast<const regionCoupledBaseFvPatch>(p)),
    method_(UNDEFINED),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr)
{}


Foam::energyRegionCoupledFvPatchScalarField::
energyRegionCoupledFvPatchScalarField
(
    const energyRegionCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<scalar>(ptf, p, iF, mapper),
    regionCoupledPatch_(refCast<const regionCoupledBaseFvPatch>(p)),
    method_(ptf.method_),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr)
{}


Foam::energyRegionCoupledFvPatchScalarField::
energyRegionCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<scalar>(p, iF, dict),
    regionCoupledPatch_(refCast<const regionCoupledBaseFvPatch>(p)),
    method_(UNDEFINED),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr)
{

    if (!isA<regionCoupledBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << regionCoupledBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
}


Foam::energyRegionCoupledFvPatchScalarField::
energyRegionCoupledFvPatchScalarField
(
    const energyRegionCoupledFvPatchScalarField& ptf
)
:
    coupledFvPatchField<scalar>(ptf),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    method_(ptf.method_),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr)
{}


Foam::energyRegionCoupledFvPatchScalarField::
energyRegionCoupledFvPatchScalarField
(
    const energyRegionCoupledFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledFvPatchField<scalar>(ptf, iF),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    method_(ptf.method_),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::energyRegionCoupledFvPatchScalarField::
snGrad() const
{
    return
        regionCoupledPatch_.patch().deltaCoeffs()
       *(*this - patchInternalField());
}


Foam::tmp<Foam::scalarField> Foam::energyRegionCoupledFvPatchScalarField::
snGrad(const scalarField&) const
{
    return snGrad();
}


void Foam::energyRegionCoupledFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    label patchi = patch().index();
    const scalarField& pp =  thermoPtr_->p().boundaryField()[patchi];

    const scalarField lWeights(weights());

    scalarField::operator=
    (
        thermoPtr_->he
        (
            pp,
            lWeights*patchInternalTemperatureField()
         + (1.0 - lWeights)*patchNeighbourTemperatureField(),
            patchi
        )
    );

    fvPatchScalarField::evaluate();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::energyRegionCoupledFvPatchScalarField::
patchNeighbourField() const
{
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

    const labelUList& nbrFaceCells = nbrPatch.faceCells();

    setMethod();

    const scalarField nbrIntT
    (
        nbrThermoPtr_->T().primitiveField(), nbrFaceCells
    );

    scalarField intNbrT
    (
        regionCoupledPatch_.regionCoupledPatch().interpolate(nbrIntT)
    );

    label patchi = patch().index();
    const scalarField& pp =  thermoPtr_->p().boundaryField()[patchi];
    tmp<scalarField> tmyHE = thermoPtr_->he(pp, intNbrT, patchi);

    return tmyHE;
}


Foam::tmp<Foam::scalarField> Foam::energyRegionCoupledFvPatchScalarField::
patchNeighbourTemperatureField() const
{
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

    const labelUList& nbrFaceCells = nbrPatch.faceCells();

    const scalarField nbrIntT
    (
        nbrThermoPtr_->T().primitiveField(), nbrFaceCells
    );

     tmp<scalarField> tintNbrT =
        regionCoupledPatch_.regionCoupledPatch().interpolate(nbrIntT);

    return tintNbrT;
}


Foam::tmp<Foam::scalarField> Foam::energyRegionCoupledFvPatchScalarField::
patchInternalTemperatureField() const
{
    const labelUList& faceCells = regionCoupledPatch_.faceCells();

    tmp<scalarField> tintT
    (
        new scalarField(thermoPtr_->T().primitiveField(), faceCells)
    );

    return tintT;
}


void Foam::energyRegionCoupledFvPatchScalarField::updateInterfaceMatrix
(
    Field<scalar>& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    setMethod();

    scalarField myHE(this->size());

    if (&psiInternal == &primitiveField())
    {
        label patchi = this->patch().index();
        const scalarField& pp =  thermoPtr_->p().boundaryField()[patchi];
        const scalarField& Tp =  thermoPtr_->T().boundaryField()[patchi];

        myHE = thermoPtr_->he(pp, Tp, patchi);
    }
    else
    {
        // NOTE: This is not correct for preconditioned solvers
        // psiInternal is not the information needed of the slave
        forAll(*this, facei)
        {
            myHE[facei] = psiInternal[regionCoupledPatch_.faceCells()[facei]];
        }
    }

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = regionCoupledPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*myHE[elemI];
    }
}


void Foam::energyRegionCoupledFvPatchScalarField::updateInterfaceMatrix
(
    Field<scalar>& result,
    const Field<scalar>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    NotImplemented;
}


void Foam::energyRegionCoupledFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        energyRegionCoupledFvPatchScalarField
    );
};


// ************************************************************************* //
