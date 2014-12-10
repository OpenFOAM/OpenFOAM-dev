/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

//- comparison operator for probes class
template<class T>
class isNotEqOp
{
public:

    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-VGREAT*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note:chould check for y != unsetVal but multiple sample cells
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
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    Field<Type> values(sample(vField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[vField.name()];

        os  << setw(w) << vField.time().timeToUserTime(vField.time().value());

        forAll(values, probeI)
        {
            os  << ' ' << setw(w) << values[probeI];
        }
        os  << endl;
    }
}


template<class Type>
void Foam::probes::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    Field<Type> values(sample(sField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[sField.name()];

        os  << setw(w) << sField.time().timeToUserTime(sField.time().value());

        forAll(values, probeI)
        {
            os  << ' ' << setw(w) << values[probeI];
        }
        os  << endl;
    }
}


template<class Type>
void Foam::probes::sampleAndWrite(const fieldGroup<Type>& fields)
{
    forAll(fields, fieldI)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fields[fieldI],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = mesh_.find(fields[fieldI]);

            if
            (
                iter != objectRegistry::end()
             && iter()->type()
             == GeometricField<Type, fvPatchField, volMesh>::typeName
            )
            {
                sampleAndWrite
                (
                    mesh_.lookupObject
                    <GeometricField<Type, fvPatchField, volMesh> >
                    (
                        fields[fieldI]
                    )
                );
            }
        }
    }
}


template<class Type>
void Foam::probes::sampleAndWriteSurfaceFields(const fieldGroup<Type>& fields)
{
    forAll(fields, fieldI)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    IOobject
                    (
                        fields[fieldI],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = mesh_.find(fields[fieldI]);

            if
            (
                iter != objectRegistry::end()
             && iter()->type()
             == GeometricField<Type, fvsPatchField, surfaceMesh>::typeName
            )
            {
                sampleAndWrite
                (
                    mesh_.lookupObject
                    <GeometricField<Type, fvsPatchField, surfaceMesh> >
                    (
                        fields[fieldI]
                    )
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::probes::sample
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    tmp<Field<Type> > tValues
    (
        new Field<Type>(this->size(), unsetVal)
    );

    Field<Type>& values = tValues();

    if (fixedLocations_)
    {
        autoPtr<interpolation<Type> > interpolator
        (
            interpolation<Type>::New(interpolationScheme_, vField)
        );

        forAll(*this, probeI)
        {
            if (elementList_[probeI] >= 0)
            {
                const vector& position = operator[](probeI);

                values[probeI] = interpolator().interpolate
                (
                    position,
                    elementList_[probeI],
                    -1
                );
            }
        }
    }
    else
    {
        forAll(*this, probeI)
        {
            if (elementList_[probeI] >= 0)
            {
                values[probeI] = vField[elementList_[probeI]];
            }
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::probes::sample(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::probes::sample
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    tmp<Field<Type> > tValues
    (
        new Field<Type>(this->size(), unsetVal)
    );

    Field<Type>& values = tValues();

    forAll(*this, probeI)
    {
        if (faceList_[probeI] >= 0)
        {
            values[probeI] = sField[faceList_[probeI]];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::probes::sampleSurfaceFields(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<GeometricField<Type, fvsPatchField, surfaceMesh> >
        (
            fieldName
        )
    );
}

// ************************************************************************* //
