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

#include "probes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOmanip.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
class isNotEqOp
{
public:

    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-vGreat*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note: should check for y != unsetVal but multiple sample cells
            // already handled in read().
        }
        else
        {
            // x is not set. y might be.
            x = y;
        }
    }
};

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::probes::sampleAndWrite
(
    const VolField<Type>& vField
)
{
    Field<Type> values(sample(vField));

    if (Pstream::master())
    {
        const unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[vField.name()];

        os  << setw(w) << vField.time().userTimeValue();

        forAll(values, probei)
        {
            OStringStream buf;
            buf << values[probei];
            os  << ' ' << setw(w) << buf.str().c_str();
        }
        os  << endl;
    }
}


template<class Type>
void Foam::probes::sampleAndWrite
(
    const SurfaceField<Type>& sField
)
{
    Field<Type> values(sample(sField));

    if (Pstream::master())
    {
        const unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[sField.name()];

        os  << sField.time().userTimeValue();

        forAll(values, probei)
        {
            OStringStream buf;
            buf << values[probei];
            os  << ' ' << setw(w) << buf.str().c_str();
        }
        os  << endl;
    }
}


template<class Type>
void Foam::probes::sampleAndWrite(const fieldGroup<Type>& fields)
{
    forAll(fields, fieldi)
    {
        objectRegistry::const_iterator iter = mesh_.find(fields[fieldi]);

        if
        (
            iter != objectRegistry::end()
         && iter()->type()
         == VolField<Type>::typeName
        )
        {
            sampleAndWrite
            (
                mesh_.lookupObject
                <VolField<Type>>
                (
                    fields[fieldi]
                )
            );
        }
    }
}


template<class Type>
void Foam::probes::sampleAndWriteSurfaceFields(const fieldGroup<Type>& fields)
{
    forAll(fields, fieldi)
    {
        objectRegistry::const_iterator iter = mesh_.find(fields[fieldi]);

        if
        (
            iter != objectRegistry::end()
         && iter()->type()
         == SurfaceField<Type>::typeName
        )
        {
            sampleAndWrite
            (
                mesh_.lookupObject
                <SurfaceField<Type>>
                (
                    fields[fieldi]
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample
(
    const VolField<Type>& vField
) const
{
    const Type unsetVal(-vGreat*pTraits<Type>::one);

    tmp<Field<Type>> tValues
    (
        new Field<Type>(this->size(), unsetVal)
    );

    Field<Type>& values = tValues.ref();

    if (fixedLocations_)
    {
        autoPtr<interpolation<Type>> interpolator
        (
            interpolation<Type>::New(interpolationScheme_, vField)
        );

        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                const vector& position = operator[](probei);

                values[probei] = interpolator().interpolate
                (
                    position,
                    elementList_[probei],
                    -1
                );
            }
        }
    }
    else
    {
        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                values[probei] = vField[elementList_[probei]];
            }
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<VolField<Type>>
        (
            fieldName
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample
(
    const SurfaceField<Type>& sField
) const
{
    const Type unsetVal(-vGreat*pTraits<Type>::one);

    tmp<Field<Type>> tValues
    (
        new Field<Type>(this->size(), unsetVal)
    );

    Field<Type>& values = tValues.ref();

    forAll(*this, probei)
    {
        if (faceList_[probei] >= 0)
        {
            values[probei] = sField[faceList_[probei]];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sampleSurfaceFields(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<SurfaceField<Type>>
        (
            fieldName
        )
    );
}

// ************************************************************************* //
