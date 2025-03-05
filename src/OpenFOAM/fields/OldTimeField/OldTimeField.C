/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "OldTimeField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
const FieldType& Foam::OldTimeField<FieldType>::field() const
{
    return static_cast<const FieldType&>(*this);
}


template<class FieldType>
FieldType& Foam::OldTimeField<FieldType>::fieldRef()
{
    return static_cast<FieldType&>(*this);
}


template<class FieldType>
void Foam::OldTimeField<FieldType>::storeOldTimesInner() const
{
    if (tfield0_.valid())
    {
        if (notNull(tfield0_()))
        {
            // Propagate to store the old-old field
            tfield0_.ref().OldTimeField<Field0Type>::storeOldTimesInner();

            // Set the old-field to this field
            tfield0_.ref() = field();
            tfield0_.ref().OldTimeField<Field0Type>::timeIndex_ = timeIndex_;

            // If we have an old-old field, then the old field is state and
            // should be written in the same way as the field
            if (tfield0_().OldTimeField<Field0Type>::tfield0_.valid())
            {
                tfield0_.ref().writeOpt() = field().writeOpt();
            }
        }
        else
        {
            // Reinstate the old-time field
            oldTime();
        }
    }
}


template<class FieldType>
void Foam::OldTimeField<FieldType>::nullOldestTimeInner()
{
    if (tfield0_.valid() && notNull(tfield0_()))
    {
        if (tfield0_().OldTimeField<Field0Type>::tfield0_.valid())
        {
            tfield0_.ref().OldTimeField<Field0Type>::nullOldestTimeInner();
        }
        else
        {
            tfield0_ = tmp<Field0Type>(NullObjectRef<FieldType>());
        }
    }
}


template<class FieldType>
template<class OldTimeBaseField>
void Foam::OldTimeField<FieldType>::setBase(const OldTimeBaseField& otbf) const
{
    if (!tfield0_.valid())
    {
        otbf.tfield0_.clear();
    }
    else
    {
        otbf.tfield0_ = tmp<typename Field0Type::Base>(tfield0_());
    }

    otbf.timeIndex_ = timeIndex_;

    otbf.setBase();
}


template<class FieldType>
void Foam::OldTimeField<FieldType>::setBase(const nil&) const
{}


template<class FieldType>
void Foam::OldTimeField<FieldType>::setBase() const
{
    setBase(OldTimeBaseFieldType<FieldType>()(*this));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class FieldType>
bool Foam::OldTimeField<FieldType>::readOldTimeIfPresent()
{
    typeIOobject<Field0Type> io0
    (
        field().name() + "_0",
        field().time().name(),
        field().db(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        field().registerObject()
    );

    if (io0.headerOk())
    {
        tfield0_ = new Field0Type(io0, field().mesh());
        setBase();

        tfield0_.ref().OldTimeField<Field0Type>::timeIndex_ = timeIndex_ - 1;
        tfield0_.ref().OldTimeField<Field0Type>::setBase();

        if (!tfield0_.ref().OldTimeField<Field0Type>::readOldTimeIfPresent())
        {
            tfield0_.ref().OldTimeField<Field0Type>::oldTime();
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class FieldType>
template<template<class> class OtherPrimitiveField>
void Foam::OldTimeField<FieldType>::copyOldTimes
(
    const IOobject& io,
    const OtherOldTime<OtherPrimitiveField>& otf
)
{
    copyOldTimes(io.name(), otf);
}


template<class FieldType>
template<template<class> class OtherPrimitiveField>
void Foam::OldTimeField<FieldType>::copyOldTimes
(
    const word& newName,
    const OtherOldTime<OtherPrimitiveField>& otf
)
{
    if (otf.tfield0_.valid() && notNull(otf.tfield0_()))
    {
        tfield0_ = new Field0Type(newName + "_0", otf.tfield0_());
        setBase();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class FieldType>
Foam::OldTimeField<FieldType>::OldTimeField(const label timeIndex)
:
    timeIndex_(timeIndex),
    tfield0_(nullptr)
{}


template<class FieldType>
Foam::OldTimeField<FieldType>::OldTimeField(const OldTimeField<FieldType>& otf)
:
    timeIndex_(otf.timeIndex_),
    tfield0_(nullptr)
{
    if (otf.tfield0_.valid() && notNull(otf.tfield0_()))
    {
        tfield0_ = new Field0Type(otf.tfield0_());
        setBase();
    }
}


template<class FieldType>
Foam::OldTimeField<FieldType>::OldTimeField(OldTimeField<FieldType>&& otf)
:
    timeIndex_(otf.timeIndex_),
    tfield0_(nullptr)
{
    if (otf.tfield0_.valid() && notNull(otf.tfield0_()))
    {
        tfield0_ = tmp<Field0Type>(move(otf.tfield0_));
        setBase();
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class FieldType>
Foam::OldTimeField<FieldType>::~OldTimeField()
{
    if (tfield0_.valid() && notNull(tfield0_()))
    {
        tfield0_.clear();
        setBase();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class FieldType>
bool Foam::OldTimeField<FieldType>::isOldTime() const
{
    return
        field().name().size() > 2
     && field().name()(field().name().size() - 2, 2) == "_0";
}


template<class FieldType>
bool Foam::OldTimeField<FieldType>::hasStoredOldTimes() const
{
    return timeIndex_ == field().time().timeIndex();
}


template<class FieldType>
void Foam::OldTimeField<FieldType>::storeOldTimes() const
{
    // Store if it has not yet been done in this time-step
    if (!hasStoredOldTimes())
    {
        // Propagate through the old-time fields
        if (!isOldTime())
        {
            storeOldTimesInner();
        }

        // Update the time index
        timeIndex_ = field().time().timeIndex();
        setBase();
    }
}


template<class FieldType>
void Foam::OldTimeField<FieldType>::clearOldTimes()
{
    if (tfield0_.valid())
    {
        tfield0_ = tmp<Field0Type>();
        setBase();
    }
}


template<class FieldType>
void Foam::OldTimeField<FieldType>::nullOldestTime()
{
    if (!isOldTime())
    {
        nullOldestTimeInner();
    }
}


template<class FieldType>
Foam::label Foam::OldTimeField<FieldType>::nOldTimes
(
    const bool includeNull
) const
{
    if (tfield0_.valid())
    {
        if (isNull(tfield0_()))
        {
            return includeNull;
        }
        else
        {
            return tfield0_->nOldTimes(includeNull) + 1;
        }
    }
    else
    {
        return 0;
    }
}


template<class FieldType>
const typename Foam::OldTimeField<FieldType>::Field0Type&
Foam::OldTimeField<FieldType>::oldTime() const
{
    if (!tfield0_.valid() || isNull(tfield0_()))
    {
        // Old-time field does not yet exist. Create it.

        // Clear the field0Ptr to ensure the old-time field constructor
        // does not construct the old-old-time field
        tfield0_.clear();
        setBase();

        // Construct a copy of the field
        tfield0_ =
            OldTimeFieldCopy<FieldType>()
            (
                IOobject
                (
                    field().name() + "_0",
                    field().time().name(),
                    field().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    field().registerObject()
                ),
                field()
            );
        setBase();
    }
    else
    {
        // Old-time field exists. Update as necessary.
        storeOldTimes();
    }

    return tfield0_();
}


template<class FieldType>
typename Foam::OldTimeField<FieldType>::Field0Type&
Foam::OldTimeField<FieldType>::oldTimeRef()
{
    static_cast<const OldTimeField<FieldType>&>(*this).oldTime();

    // Note: Not tfield0_.ref(), because this might be a base field storing a
    // reference only. It's valid to un-const this reference because the
    // derived field is guaranteed to store the non-const pointer.
    return const_cast<Field0Type&>(tfield0_());
}


template<class FieldType>
const typename Foam::OldTimeField<FieldType>::Field0Type&
Foam::OldTimeField<FieldType>::oldTime
(
    const label n
) const
{
    return n == 0 ? field() : oldTime().oldTime(n - 1);
}


template<class FieldType>
typename Foam::OldTimeField<FieldType>::Field0Type&
Foam::OldTimeField<FieldType>::oldTimeRef
(
    const label n
)
{
    return n == 0 ? fieldRef() : oldTimeRef().oldTimeRef(n - 1);
}


// ************************************************************************* //
