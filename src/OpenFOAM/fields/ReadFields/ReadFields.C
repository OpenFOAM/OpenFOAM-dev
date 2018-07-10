/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ReadFields.H"
#include "IOobjectList.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField, class Mesh>
Foam::wordList Foam::ReadFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar
)
{
    // Search list of objects for wanted type
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    wordList masterNames(fieldObjects.names());

    if (syncPar && Pstream::parRun())
    {
        // Check that I have the same fields as the master
        const wordList localNames(masterNames);
        Pstream::scatter(masterNames);

        HashSet<word> localNamesSet(localNames);

        forAll(masterNames, i)
        {
            const word& masterFld = masterNames[i];

            HashSet<word>::iterator iter = localNamesSet.find(masterFld);

            if (iter == localNamesSet.end())
            {
                FatalErrorInFunction
                    << "Fields not synchronised across processors." << endl
                    << "Master has fields " << masterNames
                    << "  processor " << Pstream::myProcNo()
                    << " has fields " << localNames << exit(FatalError);
            }
            else
            {
                localNamesSet.erase(iter);
            }
        }

        forAllConstIter(HashSet<word>, localNamesSet, iter)
        {
            FatalErrorInFunction
                << "Fields not synchronised across processors." << endl
                << "Master has fields " << masterNames
                << "  processor " << Pstream::myProcNo()
                << " has fields " << localNames << exit(FatalError);
        }
    }


    fields.setSize(masterNames.size());

    // Make sure to read in masterNames order.

    forAll(masterNames, i)
    {
        Info<< "Reading " << GeoField::typeName << ' ' << masterNames[i]
            << endl;

        const IOobject& io = *fieldObjects[masterNames[i]];

        fields.set
        (
            i,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    io.registerObject()
                ),
                mesh
            )
        );
    }
    return masterNames;
}


template<class GeoField>
void Foam::ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    objectRegistry& fieldsCache
)
{
    // Collect all times that are no longer used
    {
        HashSet<word> usedTimes(timeNames);

        DynamicList<word> unusedTimes(fieldsCache.size());

        forAllIter(objectRegistry, fieldsCache, timeIter)
        {
            const word& tm = timeIter.key();
            if (!usedTimes.found(tm))
            {
                unusedTimes.append(tm);
            }
        }

        forAll(unusedTimes, i)
        {
            objectRegistry& timeCache = const_cast<objectRegistry&>
            (
                fieldsCache.lookupObject<objectRegistry>(unusedTimes[i])
            );
            fieldsCache.checkOut(timeCache);
        }
    }


    // Load any new fields
    forAll(timeNames, i)
    {
        const word& tm = timeNames[i];

        // Create if not found
        if (!fieldsCache.found(tm))
        {
            // Create objectRegistry if not found
            objectRegistry* timeCachePtr = new objectRegistry
            (
                IOobject
                (
                    tm,
                    tm,
                    fieldsCache,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );
            timeCachePtr->store();
        }

        // Obtain cache for current time
        const objectRegistry& timeCache =
            fieldsCache.lookupObject<objectRegistry>
            (
                tm
            );

        // Store field if not found
        if (!timeCache.found(fieldName))
        {
            GeoField loadedFld
            (
                IOobject
                (
                    fieldName,
                    tm,
                    mesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            );

            // Transfer to timeCache (new objectRegistry and store flag)
            GeoField* fldPtr = new GeoField
            (
                IOobject
                (
                    fieldName,
                    tm,
                    timeCache,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                loadedFld
            );
            fldPtr->store();
        }
    }
}


template<class GeoField>
void Foam::ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    const word& registryName
)
{
    ReadFields<GeoField>
    (
        fieldName,
        mesh,
        timeNames,
        const_cast<objectRegistry&>
        (
            mesh.thisDb().subRegistry(registryName, true)
        )
    );
}


template<class GeoFieldType>
void Foam::readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    LIFOStack<regIOobject*>& storedObjects
)
{
    IOobjectList fields(objects.lookupClass(GeoFieldType::typeName));
    if (!fields.size()) return;

    bool firstField = true;

    forAllConstIter(IOobjectList, fields, fieldIter)
    {
        const IOobject& io = *fieldIter();
        const word& fieldName = io.name();

        if (selectedFields.found(fieldName))
        {
            if (firstField)
            {
                Info<< "    " << GeoFieldType::typeName << "s:";
                firstField = false;
            }

            Info<< " " << fieldName;

            GeoFieldType* fieldPtr = new GeoFieldType
            (
                IOobject
                (
                    fieldName,
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            fieldPtr->store();
            storedObjects.push(fieldPtr);
        }
    }

    if (!firstField)
    {
        Info<< endl;
    }
}


template<class UniformFieldType>
void Foam::readUniformFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    LIFOStack<regIOobject*>& storedObjects,
    const bool syncPar
)
{
    // Search list of objects for wanted type
    IOobjectList fields(objects.lookupClass(UniformFieldType::typeName));
    if (!fields.size()) return;

    wordList masterNames(fields.names());

    if (syncPar && Pstream::parRun())
    {
        // Check that I have the same fields as the master
        const wordList localNames(masterNames);
        Pstream::scatter(masterNames);

        HashSet<word> localNamesSet(localNames);

        forAll(masterNames, i)
        {
            const word& masterFld = masterNames[i];

            HashSet<word>::iterator iter = localNamesSet.find(masterFld);

            if (iter == localNamesSet.end())
            {
                FatalErrorInFunction
                    << "Fields not synchronised across processors." << endl
                    << "Master has fields " << masterNames
                    << "  processor " << Pstream::myProcNo()
                    << " has fields " << localNames << exit(FatalError);
            }
            else
            {
                localNamesSet.erase(iter);
            }
        }

        forAllConstIter(HashSet<word>, localNamesSet, iter)
        {
            FatalErrorInFunction
                << "Fields not synchronised across processors." << endl
                << "Master has fields " << masterNames
                << "  processor " << Pstream::myProcNo()
                << " has fields " << localNames << exit(FatalError);
        }
    }

    bool firstField = true;

    forAll(masterNames, i)
    {
        const IOobject& io = *fields[masterNames[i]];
        const word& fieldName = io.name();

        if (selectedFields.found(fieldName))
        {
            if (firstField)
            {
                Info<< "    " << UniformFieldType::typeName << "s:";
                firstField = false;
            }

            Info<< " " << fieldName;

            UniformFieldType* fieldPtr = new UniformFieldType
            (
                IOobject
                (
                    fieldName,
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            fieldPtr->store();
            storedObjects.push(fieldPtr);
        }
    }

    Info<< endl;
}


// ************************************************************************* //
