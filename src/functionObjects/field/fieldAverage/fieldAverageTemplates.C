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

#include "fieldAverageItem.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldAverage::addMeanFieldType(const label fieldi)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    Log << "    Reading/initialising field " << meanFieldName << endl;

    if (obr_.foundObject<Type>(meanFieldName))
    {}
    else if (obr_.found(meanFieldName))
    {
        Log << "    Cannot allocate average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldi].mean() = false;
    }
    else
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        // Store on registry
        obr_.store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr_.time().timeName(obr_.time().startTime().value()),
                    obr_,
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::addMeanField(const label fieldi)
{
    if (faItems_[fieldi].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh>
            VolFieldType;

        typedef typename VolFieldType::Internal InternalType;

        typedef GeometricField<Type, fvsPatchField, surfaceMesh>
            SurfaceFieldType;

        const word& fieldName = faItems_[fieldi].fieldName();

        if (obr_.foundObject<VolFieldType>(fieldName))
        {
            addMeanFieldType<VolFieldType>(fieldi);
        }
        else if (obr_.foundObject<InternalType>(fieldName))
        {
            addMeanFieldType<InternalType>(fieldi);
        }
        else if (obr_.foundObject<SurfaceFieldType>(fieldName))
        {
            addMeanFieldType<SurfaceFieldType>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addPrime2MeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    Log << "    Reading/initialising field " << prime2MeanFieldName << nl;

    if (obr_.foundObject<Type2>(prime2MeanFieldName))
    {}
    else if (obr_.found(prime2MeanFieldName))
    {
        Log << "    Cannot allocate average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField = obr_.lookupObject<Type1>(meanFieldName);

        // Store on registry
        obr_.store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr_.time().timeName(obr_.time().startTime().value()),
                    obr_,
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addPrime2MeanField(const label fieldi)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef typename VolFieldType1::Internal InternalType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef typename VolFieldType2::Internal InternalType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (obr_.foundObject<VolFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<VolFieldType1, VolFieldType2>(fieldi);
        }
        else if (obr_.foundObject<InternalType1>(fieldName))
        {
            addPrime2MeanFieldType<InternalType1, InternalType2>(fieldi);
        }
        else if (obr_.foundObject<SurfaceFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>
            (
                fieldi
            );
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        Type& meanField =
            obr_.lookupObjectRef<Type>(faItems_[fieldi].meanFieldName());

        scalar dt = obr_.time().deltaTValue();
        scalar Dt = totalTime_[fieldi];

        if (iterBase())
        {
            dt = 1;
            Dt = scalar(totalIter_[fieldi]);
        }

        scalar beta = dt/Dt;

        if (window() > 0)
        {
            const scalar w = window();

            if (Dt - dt >= w)
            {
                beta = dt/w;
            }
        }

        meanField = (1 - beta)*meanField + beta*baseField;
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal InternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr_.foundObject<VolFieldType>(fieldName))
            {
                calculateMeanFieldType<VolFieldType>(fieldi);
            }
            else if (obr_.foundObject<InternalType>(fieldName))
            {
                calculateMeanFieldType<InternalType>(fieldi);
            }
            else if (obr_.foundObject<SurfaceFieldType>(fieldName))
            {
                calculateMeanFieldType<SurfaceFieldType>(fieldi);
            }
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
    const Type1& meanField =
        obr_.lookupObject<Type1>(faItems_[fieldi].meanFieldName());

    Type2& prime2MeanField =
        obr_.lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

    scalar dt = obr_.time().deltaTValue();
    scalar Dt = totalTime_[fieldi];

    if (iterBase())
    {
        dt = 1;
        Dt = scalar(totalIter_[fieldi]);
    }

    scalar beta = dt/Dt;

    if (window() > 0)
    {
        const scalar w = window();

        if (Dt - dt >= w)
        {
            beta = dt/w;
        }
    }

    prime2MeanField =
        (1 - beta)*prime2MeanField
      + beta*sqr(baseField)
      - sqr(meanField);
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFields() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef typename VolFieldType1::Internal InternalType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef typename VolFieldType2::Internal InternalType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr_.foundObject<VolFieldType1>(fieldName))
            {
                calculatePrime2MeanFieldType<VolFieldType1, VolFieldType2>
                (
                    fieldi
                );
            }
            else if (obr_.foundObject<InternalType1>(fieldName))
            {
                calculatePrime2MeanFieldType<InternalType1, InternalType2>
                (
                    fieldi
                );
            }
            else if (obr_.foundObject<SurfaceFieldType1>(fieldName))
            {
                calculatePrime2MeanFieldType
                <
                    SurfaceFieldType1,
                    SurfaceFieldType2
                >(fieldi);
            }
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2MeanType
(
    const label fieldi
) const
{
    const Type1& meanField =
        obr_.lookupObject<Type1>(faItems_[fieldi].meanFieldName());

    Type2& prime2MeanField =
        obr_.lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

    prime2MeanField += sqr(meanField);
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2Mean() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef typename VolFieldType1::Internal InternalType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef typename VolFieldType2::Internal InternalType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr_.foundObject<VolFieldType1>(fieldName))
            {
                addMeanSqrToPrime2MeanType<VolFieldType1, VolFieldType2>
                (
                    fieldi
                );
            }
            else if (obr_.foundObject<InternalType1>(fieldName))
            {
                addMeanSqrToPrime2MeanType<InternalType1, InternalType2>
                (
                    fieldi
                );
            }
            else if (obr_.foundObject<SurfaceFieldType1>(fieldName))
            {
                addMeanSqrToPrime2MeanType<SurfaceFieldType1, SurfaceFieldType2>
                (
                    fieldi
                );
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFieldType
(
    const word& fieldName
) const
{
    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& f = obr_.lookupObject<Type>(fieldName);
        f.write();
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal InternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& fieldName = faItems_[fieldi].meanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<InternalType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
        }
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].prime2MeanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<InternalType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
        }
    }
}


// ************************************************************************* //
