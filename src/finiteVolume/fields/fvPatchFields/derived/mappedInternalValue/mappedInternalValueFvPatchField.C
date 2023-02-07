/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "mappedInternalValueFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::mappedInternalPatchBase&
Foam::mappedInternalValueFvPatchField<Type>::mapper() const
{
    if (mapperPtr_.valid())
    {
        return mapperPtr_();
    }

    if (isA<mappedInternalPatchBase>(this->patch().patch()))
    {
        return refCast<const mappedInternalPatchBase>(this->patch().patch());
    }

    FatalErrorInFunction
        << "Field " << this->internalField().name() << " on patch "
        << this->patch().name() << " in file "
        << this->internalField().objectPath()
        << " has neither a mapper specified nor is the patch of "
        << mappedInternalPatchBase::typeName << " type"
        << exit(FatalError);

    return NullObjectRef<mappedInternalPatchBase>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedInternalValueFvPatchField<Type>::
mappedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldName_(iF.name()),
    setAverage_(false),
    average_(Zero),
    interpolationScheme_(interpolationCell<Type>::typeName),
    mapperPtr_(nullptr)
{}


template<class Type>
Foam::mappedInternalValueFvPatchField<Type>::
mappedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    fieldName_(dict.lookupOrDefault<word>("field", iF.name())),
    setAverage_
    (
        dict.lookupOrDefault<bool>("setAverage", dict.found("average"))
    ),
    average_(setAverage_ ? dict.lookup<Type>("average") : Zero),
    interpolationScheme_(dict.lookup<word>("interpolationScheme")),
    mapperPtr_
    (
        mappedInternalPatchBase::specified(dict)
      ? new mappedInternalPatchBase(p.patch(), dict)
      : nullptr
    )
{}


template<class Type>
Foam::mappedInternalValueFvPatchField<Type>::
mappedInternalValueFvPatchField
(
    const mappedInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_),
    interpolationScheme_(ptf.interpolationScheme_),
    mapperPtr_
    (
        ptf.mapperPtr_.valid()
      ? new mappedInternalPatchBase(p.patch(), ptf.mapperPtr_())
      : nullptr
    )
{}


template<class Type>
Foam::mappedInternalValueFvPatchField<Type>::
mappedInternalValueFvPatchField
(
    const mappedInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldName_(ptf.fieldName_),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_),
    interpolationScheme_(ptf.interpolationScheme_),
    mapperPtr_
    (
        ptf.mapperPtr_.valid()
      ? new mappedInternalPatchBase(ptf.patch().patch(), ptf.mapperPtr_())
      : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedInternalValueFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::map(ptf, mapper);

    if (mapperPtr_.valid())
    {
        mapperPtr_->clearOut();
    }
}


template<class Type>
void Foam::mappedInternalValueFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    fixedValueFvPatchField<Type>::reset(ptf);

    if (mapperPtr_.valid())
    {
        mapperPtr_->clearOut();
    }
}


template<class Type>
void Foam::mappedInternalValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const fvMesh& nbrMesh = refCast<const fvMesh>(this->mapper().nbrMesh());

    const VolField<Type>& nbrField =
        this->mapper().sameRegion()
     && this->fieldName_ == this->internalField().name()
      ? refCast<const VolField<Type>>(this->internalField())
      : nbrMesh.template lookupObject<VolField<Type>>(this->fieldName_);

    // Construct mapped values
    Field<Type> sampleValues;

    if (interpolationScheme_ != interpolationCell<Type>::typeName)
    {
        // Create an interpolation
        autoPtr<interpolation<Type>> interpolatorPtr
        (
            interpolation<Type>::New
            (
                interpolationScheme_,
                nbrField
            )
        );
        const interpolation<Type>& interpolator = interpolatorPtr();

        // Cells on which samples are generated
        const labelList& sampleCells = mapper().cellIndices();

        // Send the patch points to the cells
        pointField samplePoints(mapper().samplePoints());
        mapper().map().reverseDistribute
        (
            sampleCells.size(),
            samplePoints
        );

        // Interpolate values
        sampleValues.resize(sampleCells.size());
        forAll(sampleCells, i)
        {
            if (sampleCells[i] != -1)
            {
                sampleValues[i] =
                    interpolator.interpolate
                    (
                        samplePoints[i],
                        sampleCells[i]
                    );
            }
        }

        // Send the values back to the patch
        mapper().map().distribute(sampleValues);
    }
    else
    {
        // No interpolation. Just sample cell values directly.
        sampleValues = mapper().distribute(nbrField);
    }

    // Set the average, if necessary
    if (setAverage_)
    {
        const Type sampleAverageValue =
            gSum(this->patch().magSf()*sampleValues)
           /gSum(this->patch().magSf());

        if (mag(sampleAverageValue)/mag(average_) > 0.5)
        {
            sampleValues *= mag(average_)/mag(sampleAverageValue);
        }
        else
        {
            sampleValues += average_ - sampleAverageValue;
        }
    }

    // Assign sampled patch values
    this->operator==(sampleValues);

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedInternalValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    writeEntryIfDifferent
    (
        os,
        "field",
        this->internalField().name(),
        fieldName_
    );

    if (setAverage_)
    {
        writeEntry(os, "average", average_);
    }

    writeEntry(os, "interpolationScheme", interpolationScheme_);

    if (mapperPtr_.valid())
    {
        mapperPtr_->write(os);
    }

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
