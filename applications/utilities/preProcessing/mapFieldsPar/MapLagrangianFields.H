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

InNamespace
    Foam

Description
    Gets the indices of (source)particles that have been appended to the
    target cloud and maps the lagrangian fields accordingly.

\*---------------------------------------------------------------------------*/

#ifndef MapLagrangianFields_H
#define MapLagrangianFields_H

#include "cloud.H"
#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"
#include "CompactIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Gets the indices of (source)particles that have been appended to the
//  target cloud and maps the lagrangian fields accordingly.
template<class Type>
void MapLagrangianFields
(
    const string& cloudName,
    const IOobjectList& objects,
    const polyMesh& meshTarget,
    const labelList& addParticles
)
{
    {
        IOobjectList fields = objects.lookupClass(IOField<Type>::typeName);

        forAllIter(IOobjectList, fields, fieldIter)
        {
            const word& fieldName = fieldIter()->name();

            Info<< "    mapping lagrangian field " << fieldName << endl;

            // Read field (does not need mesh)
            IOField<Type> fieldSource(*fieldIter());

            // Map
            IOField<Type> fieldTarget
            (
                IOobject
                (
                    fieldName,
                    meshTarget.time().timeName(),
                    cloud::prefix/cloudName,
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                addParticles.size()
            );

            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }

            // Write field
            fieldTarget.write();
        }
    }

    {
        IOobjectList fieldFields =
            objects.lookupClass(IOField<Field<Type>>::typeName);

        forAllIter(IOobjectList, fieldFields, fieldIter)
        {
            const word& fieldName = fieldIter()->name();

            Info<< "    mapping lagrangian fieldField " << fieldName << endl;

            // Read field (does not need mesh)
            IOField<Field<Type>> fieldSource(*fieldIter());

            // Map - use CompactIOField to automatically write in
            // compact form for binary format.
            CompactIOField<Field<Type>, Type> fieldTarget
            (
                IOobject
                (
                    fieldName,
                    meshTarget.time().timeName(),
                    cloud::prefix/cloudName,
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                addParticles.size()
            );

            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }

            // Write field
            fieldTarget.write();
        }
    }

    {
        IOobjectList fieldFields =
            objects.lookupClass(CompactIOField<Field<Type>, Type>::typeName);

        forAllIter(IOobjectList, fieldFields, fieldIter)
        {
            Info<< "    mapping lagrangian fieldField "
                << fieldIter()->name() << endl;

            // Read field (does not need mesh)
            CompactIOField<Field<Type>, Type> fieldSource(*fieldIter());

            // Map
            CompactIOField<Field<Type>, Type> fieldTarget
            (
                IOobject
                (
                    fieldIter()->name(),
                    meshTarget.time().timeName(),
                    cloud::prefix/cloudName,
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                addParticles.size()
            );

            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }

            // Write field
            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
