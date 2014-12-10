/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fieldAverage::addMeanFieldType(const label fieldI)
{
    faItems_[fieldI].active() = true;

    const word& fieldName = faItems_[fieldI].fieldName();
    const word& meanFieldName = faItems_[fieldI].meanFieldName();

    Info<< "    Reading/initialising field " << meanFieldName << endl;

    if (obr_.foundObject<Type>(meanFieldName))
    {
       // do nothing
    }
    else if (obr_.found(meanFieldName))
    {
        Info<< "    Cannot allocate average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldI].mean() = false;
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
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::fieldAverage::addMeanField(const label fieldI)
{
    if (faItems_[fieldI].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
        typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfFieldType;

        const word& fieldName = faItems_[fieldI].fieldName();

        if (obr_.foundObject<volFieldType>(fieldName))
        {
            addMeanFieldType<volFieldType>(fieldI);
        }
        else if (obr_.foundObject<surfFieldType>(fieldName))
        {
            addMeanFieldType<surfFieldType>(fieldI);
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addPrime2MeanFieldType(const label fieldI)
{
    const word& fieldName = faItems_[fieldI].fieldName();
    const word& meanFieldName = faItems_[fieldI].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldI].prime2MeanFieldName();

    Info<< "    Reading/initialising field " << prime2MeanFieldName << nl;

    if (obr_.foundObject<Type2>(prime2MeanFieldName))
    {
        // do nothing
    }
    else if (obr_.found(prime2MeanFieldName))
    {
        Info<< "    Cannot allocate average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldI].prime2Mean() = false;
    }
    else
    {
        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField = obr_.lookupObject<Type1>(meanFieldName);

        obr_.store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr_.time().timeName(obr_.time().startTime().value()),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addPrime2MeanField(const label fieldI)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> volFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> surfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> volFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> surfFieldType2;

    if (faItems_[fieldI].prime2Mean())
    {
        const word& fieldName = faItems_[fieldI].fieldName();

        if (!faItems_[fieldI].mean())
        {
            FatalErrorIn
            (
                "void Foam::fieldAverage::addPrime2MeanField(const label) const"
            )
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (obr_.foundObject<volFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<volFieldType1, volFieldType2>(fieldI);
        }
        else if (obr_.foundObject<surfFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<surfFieldType1, surfFieldType2>(fieldI);
        }
    }
}


template<class Type>
void Foam::fieldAverage::calculateMeanFieldType(const label fieldI) const
{
    const word& fieldName = faItems_[fieldI].fieldName();

    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        Type& meanField = const_cast<Type&>
        (
            obr_.lookupObject<Type>(faItems_[fieldI].meanFieldName())
        );

        scalar dt = obr_.time().deltaTValue();
        scalar Dt = totalTime_[fieldI];

        if (faItems_[fieldI].iterBase())
        {
            dt = 1.0;
            Dt = scalar(totalIter_[fieldI]);
        }

        scalar alpha = (Dt - dt)/Dt;
        scalar beta = dt/Dt;

        if (faItems_[fieldI].window() > 0)
        {
            const scalar w = faItems_[fieldI].window();

            if (Dt - dt >= w)
            {
                alpha = (w - dt)/w;
                beta = dt/w;
            }
        }

        meanField = alpha*meanField + beta*baseField;
    }
}


template<class Type>
void Foam::fieldAverage::calculateMeanFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfFieldType;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            calculateMeanFieldType<volFieldType>(i);
            calculateMeanFieldType<surfFieldType>(i);
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::calculatePrime2MeanFieldType(const label fieldI) const
{
    const word& fieldName = faItems_[fieldI].fieldName();

    if (obr_.foundObject<Type1>(fieldName))
    {
        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField =
            obr_.lookupObject<Type1>(faItems_[fieldI].meanFieldName());

        Type2& prime2MeanField = const_cast<Type2&>
        (
            obr_.lookupObject<Type2>(faItems_[fieldI].prime2MeanFieldName())
        );

        scalar dt = obr_.time().deltaTValue();
        scalar Dt = totalTime_[fieldI];

        if (faItems_[fieldI].iterBase())
        {
            dt = 1.0;
            Dt = scalar(totalIter_[fieldI]);
        }

        scalar alpha = (Dt - dt)/Dt;
        scalar beta = dt/Dt;

        if (faItems_[fieldI].window() > 0)
        {
            const scalar w = faItems_[fieldI].window();

            if (Dt - dt >= w)
            {
                alpha = (w - dt)/w;
                beta = dt/w;
            }
        }

        prime2MeanField =
            alpha*prime2MeanField
          + beta*sqr(baseField)
          - sqr(meanField);
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::calculatePrime2MeanFields() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> volFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> surfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> volFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> surfFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            calculatePrime2MeanFieldType<volFieldType1, volFieldType2>(i);
            calculatePrime2MeanFieldType<surfFieldType1, surfFieldType2>(i);
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addMeanSqrToPrime2MeanType(const label fieldI) const
{
    const word& fieldName = faItems_[fieldI].fieldName();

    if (obr_.foundObject<Type1>(fieldName))
    {
        const Type1& meanField =
            obr_.lookupObject<Type1>(faItems_[fieldI].meanFieldName());

        Type2& prime2MeanField = const_cast<Type2&>
        (
            obr_.lookupObject<Type2>(faItems_[fieldI].prime2MeanFieldName())
        );

        prime2MeanField += sqr(meanField);
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addMeanSqrToPrime2Mean() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> volFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> surfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> volFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> surfFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            addMeanSqrToPrime2MeanType<volFieldType1, volFieldType2>(i);
            addMeanSqrToPrime2MeanType<surfFieldType1, surfFieldType2>(i);
        }
    }
}


template<class Type>
void Foam::fieldAverage::writeFieldType(const word& fieldName) const
{
    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& f = obr_.lookupObject<Type>(fieldName);
        f.write();
    }
}


template<class Type>
void Foam::fieldAverage::writeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfFieldType;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            const word& fieldName = faItems_[i].meanFieldName();
            writeFieldType<volFieldType>(fieldName);
            writeFieldType<surfFieldType>(fieldName);
        }
        if (faItems_[i].prime2Mean())
        {
            const word& fieldName = faItems_[i].prime2MeanFieldName();
            writeFieldType<volFieldType>(fieldName);
            writeFieldType<surfFieldType>(fieldName);
        }
    }
}


// ************************************************************************* //
