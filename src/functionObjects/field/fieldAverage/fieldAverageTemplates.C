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

#include "fieldAverage.H"
#include "fieldAverageItem.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldAverage::readMeanFieldType(const label fieldi)
{
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    IOobject meanFieldIo
    (
        meanFieldName,
        obr_.time().name(),
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if
    (
        meanFieldIo.headerOk()
     && meanFieldIo.headerClassName() == Type::typeName
    )
    {
        if (obr_.found(meanFieldName))
        {
            Log << "    Cannot read average field " << meanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging for field." << endl;

            faItems_[fieldi].mean() = false;
        }
        else
        {
            Log << "    Reading field " << meanFieldName << endl;

            obr_.store(new Type(meanFieldIo, mesh_));
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::readMeanField(const label fieldi)
{
    if (faItems_[fieldi].mean())
    {
        readMeanFieldType<VolField<Type>>(fieldi);
        readMeanFieldType<VolInternalField<Type>>(fieldi);
        readMeanFieldType<SurfaceField<Type>>(fieldi);
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::initialiseMeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    if (obr_.foundObject<Type>(meanFieldName))
    {
        // Do nothing ...
    }
    else if (obr_.found(meanFieldName))
    {
        Log << "    Cannot initialise average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldi].mean() = false;
    }
    else
    {
        Log << "    Initialising field " << meanFieldName << endl;

        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        obr_.store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr_.time().name(),
                    obr_
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::initialiseMeanField
(
    const label fieldi
)
{
    if (faItems_[fieldi].mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (obr_.foundObject<VolField<Type>>(fieldName))
        {
            initialiseMeanFieldType<VolField<Type>>(fieldi);
        }
        else if (obr_.foundObject<VolInternalField<Type>>(fieldName))
        {
            initialiseMeanFieldType<VolInternalField<Type>>(fieldi);
        }
        else if (obr_.foundObject<SurfaceField<Type>>(fieldName))
        {
            initialiseMeanFieldType<SurfaceField<Type>>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::readPrime2MeanFieldType
(
    const label fieldi
)
{
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    IOobject prime2MeanFieldIo
    (
        prime2MeanFieldName,
        obr_.time().name(),
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if
    (
        prime2MeanFieldIo.headerOk()
     && prime2MeanFieldIo.headerClassName() == Type2::typeName
    )
    {
        if (obr_.found(prime2MeanFieldName))
        {
            Log << "    Cannot read average field " << prime2MeanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging for field." << endl;

            faItems_[fieldi].mean() = false;
        }
        else
        {
            Log << "    Reading field " << prime2MeanFieldName << endl;

            obr_.store(new Type2(prime2MeanFieldIo, mesh_));
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::readPrime2MeanField
(
    const label fieldi
)
{
    if (faItems_[fieldi].prime2Mean())
    {
        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << faItems_[fieldi].fieldName() << nl << exit(FatalError);
        }

        readPrime2MeanFieldType
        <VolField<Type1>, VolField<Type2>>(fieldi);
        readPrime2MeanFieldType
        <VolInternalField<Type1>, VolInternalField<Type2>>(fieldi);
        readPrime2MeanFieldType
        <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    if (obr_.foundObject<Type2>(prime2MeanFieldName))
    {
        // Do nothing ...
    }
    else if (obr_.found(prime2MeanFieldName))
    {
        Log << "    Cannot initialise average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        Log << "    Initialising field " << prime2MeanFieldName << nl;

        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField = obr_.lookupObject<Type1>(meanFieldName);

        obr_.store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr_.time().name(),
                    obr_
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanField
(
    const label fieldi
)
{
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

        if (obr_.foundObject<VolField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldType
            <VolField<Type1>, VolField<Type2>>(fieldi);
        }
        else if (obr_.foundObject<VolInternalField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldType
            <VolInternalField<Type1>, VolInternalField<Type2>>(fieldi);
        }
        else if (obr_.foundObject<SurfaceField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldType
            <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
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


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFields() const
{
    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr_.foundObject<VolField<Type>>(fieldName))
            {
                calculateMeanFieldType<VolField<Type>>(fieldi);
            }
            else if (obr_.foundObject<VolInternalField<Type>>(fieldName))
            {
                calculateMeanFieldType<VolInternalField<Type>>(fieldi);
            }
            else if (obr_.foundObject<SurfaceField<Type>>(fieldName))
            {
                calculateMeanFieldType<SurfaceField<Type>>(fieldi);
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
    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr_.foundObject<VolField<Type1>>(fieldName))
            {
                calculatePrime2MeanFieldType
                <VolField<Type1>, VolField<Type2>>(fieldi);
            }
            else if (obr_.foundObject<VolInternalField<Type1>>(fieldName))
            {
                calculatePrime2MeanFieldType
                <VolInternalField<Type1>, VolInternalField<Type2>>(fieldi);
            }
            else if (obr_.foundObject<SurfaceField<Type1>>(fieldName))
            {
                calculatePrime2MeanFieldType
                <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
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
    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr_.foundObject<VolField<Type1>>(fieldName))
            {
                addMeanSqrToPrime2MeanType
                <VolField<Type1>, VolField<Type2>>(fieldi);
            }
            else if (obr_.foundObject<VolInternalField<Type1>>(fieldName))
            {
                addMeanSqrToPrime2MeanType
                <VolInternalField<Type1>, VolInternalField<Type2>>(fieldi);
            }
            else if (obr_.foundObject<SurfaceField<Type1>>(fieldName))
            {
                addMeanSqrToPrime2MeanType
                <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
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
    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& fieldName = faItems_[fieldi].meanFieldName();
            writeFieldType<VolField<Type>>(fieldName);
            writeFieldType<VolInternalField<Type>>(fieldName);
            writeFieldType<SurfaceField<Type>>(fieldName);
        }
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].prime2MeanFieldName();
            writeFieldType<VolField<Type>>(fieldName);
            writeFieldType<VolInternalField<Type>>(fieldName);
            writeFieldType<SurfaceField<Type>>(fieldName);
        }
    }
}


// ************************************************************************* //
