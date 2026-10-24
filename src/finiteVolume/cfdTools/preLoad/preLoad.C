/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "preLoad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(preLoad, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::preLoad::preLoad(const fvMesh& mesh)
:
    timeIOdictionary
    (
        IOobject
        (
            typeName,
            mesh.time().name(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        typeName
    )
{
    if (found("fields"))
    {
        const List<word> fields(lookup("fields"));

        forAll(fields, i)
        {
            IOobject fieldIo
            (
                fields[i],
                mesh.time().name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldIo.headerOk())
            {
                if (debug)
                {
                    Info << "preLoad: Reading field " << fields[i] << endl;
                }

                #define readField(Type, GeoField)                              \
                    if (fieldIo.headerClassName() == GeoField<Type>::typeName) \
                    {                                                          \
                        regIOobject::store(new GeoField<Type>(fieldIo, mesh)); \
                        continue;                                              \
                    }
                FOR_ALL_FIELD_TYPES(readField, VolField)
                FOR_ALL_FIELD_TYPES(readField, VolInternalField);
                FOR_ALL_FIELD_TYPES(readField, SurfaceField);
            }
        }
    }
}


Foam::preLoad& Foam::preLoad::New(const fvMesh& mesh)
{
    if (mesh.db().objectRegistry::template foundObject<preLoad>(typeName))
    {
        return mesh.db().objectRegistry::template lookupObjectRef<preLoad>
        (
            typeName
        );
    }
    else
    {
        if (debug)
        {
            Info<< "preLoad::New(" << fvMesh::typeName
                << "&) : constructing " << typeName
                << " for region " << mesh.name() << endl;
        }

        preLoad* objectPtr = new preLoad(mesh);

        return regIOobject::store(objectPtr);
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::preLoad::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool updateInstance
) const
{
    if (!fields.empty())
    {
        // Add/update the fields entry in the preLoad dictionary
        const_cast<preLoad&>(*this).set("fields", fields);

        // Write the preLoad timeIOdictionary into <time>/uniform
        return timeIOdictionary::writeObject(fmt, ver, cmp, updateInstance);
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
