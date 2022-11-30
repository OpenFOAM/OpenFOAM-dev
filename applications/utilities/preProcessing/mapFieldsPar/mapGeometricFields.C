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

#include "mapGeometricFields.H"
#include "fvMeshToFvMesh.H"
#include "surfaceMesh.H"
#include "pointMesh.H"
#include "IOobjectList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void evaluateConstraintTypes(GeometricField<Type, fvPatchField, volMesh>& fld)
{
    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& fldBf = fld.boundaryFieldRef();

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(Pstream::defaultCommsType);
            }
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(Pstream::defaultCommsType);
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (patchSchedule[patchEvali].init)
                {
                    tgtField.initEvaluate(Pstream::commsTypes::scheduled);
                }
                else
                {
                    tgtField.evaluate(Pstream::commsTypes::scheduled);
                }
            }
        }
    }
}


template<class Type>
void mapVolTypeFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    const fvMeshToFvMesh& interp
)
{
    const fvMesh& srcMesh = static_cast<const fvMesh&>(interp.srcMesh());
    const fvMesh& tgtMesh = static_cast<const fvMesh&>(interp.tgtMesh());

    IOobjectList fields = objects.lookupClass(VolField<Type>::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            const VolField<Type> fieldSource(*fieldIter(), srcMesh);

            typeIOobject<VolField<Type>> targetIO
            (
                fieldName,
                tgtMesh.time().name(),
                tgtMesh,
                IOobject::MUST_READ
            );

            if (targetIO.headerOk())
            {
                Info<< "    interpolating onto existing field "
                    << fieldName << endl;
                VolField<Type> fieldTarget(targetIO, tgtMesh);

                interp.mapSrcToTgt(fieldSource, fieldTarget);

                evaluateConstraintTypes(fieldTarget);

                fieldTarget.write();
            }
            else
            {
                Info<< "    creating new field "
                    << fieldName << endl;

                targetIO.readOpt() = IOobject::NO_READ;

                tmp<VolField<Type>> tfieldTarget
                (
                    interp.mapSrcToTgt(fieldSource)
                );

                VolField<Type> fieldTarget(targetIO, tfieldTarget);

                evaluateConstraintTypes(fieldTarget);

                fieldTarget.write();
            }
        }
    }
}


template<class Type, template<class> class GeoField>
void unMappedTypeFields(const IOobjectList& objects)
{
    IOobjectList fields = objects.lookupClass(GeoField<Type>::typeName);

    forAllConstIter(IOobjectList, fields, fieldIter)
    {
        mvBak(fieldIter()->objectPath(false), "unmapped");
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::mapGeometricFields
(
    const fvMeshToFvMesh& interp,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    const polyMesh& srcMesh = interp.srcMesh();
    const polyMesh& tgtMesh = interp.tgtMesh();

    {
        // Search for list of source objects for this time
        IOobjectList objects(srcMesh, srcMesh.time().name());

        // Map the fields
        #define MapVolTypeFields(Type, nullArg)                                \
            mapVolTypeFields<Type>                                             \
            (                                                                  \
                objects,                                                       \
                selectedFields,                                                \
                interp                                                         \
            );
        FOR_ALL_FIELD_TYPES(MapVolTypeFields);
        #undef MapVolTypeFields
    }

    {
        // Search for list of target objects for this time
        IOobjectList objects(tgtMesh, tgtMesh.time().name());

        // Mark surface and point fields as unmapped
        #define UnMappedTypeFields(Type, GeoField)                             \
            unMappedTypeFields<Type, GeoField>(objects);
        FOR_ALL_FIELD_TYPES(UnMappedTypeFields, SurfaceField);
        FOR_ALL_FIELD_TYPES(UnMappedTypeFields, PointField);
        #undef UnMappedTypeFields
    }
}


// ************************************************************************* //
