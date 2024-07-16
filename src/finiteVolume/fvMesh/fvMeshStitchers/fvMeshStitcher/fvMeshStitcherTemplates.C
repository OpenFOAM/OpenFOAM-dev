/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

Description
    Perform mapping of finite volume fields required by stitching.

\*---------------------------------------------------------------------------*/

#include "fvMeshStitcher.H"
#include "conformedFvPatchField.H"
#include "conformedFvsPatchField.H"
#include "nonConformalErrorFvPatch.H"
#include "setSizeFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class GeoField>
void Foam::fvMeshStitcher::resizePatchFields()
{
    UPtrList<GeoField<Type>> fields(mesh_.fields<GeoField<Type>>());

    forAll(fields, i)
    {
        forAll(mesh_.boundary(), patchi)
        {
            typename GeoField<Type>::Patch& pf =
                boundaryFieldRefNoUpdate(fields[i])[patchi];

            if (isA<nonConformalFvPatch>(pf.patch()))
            {
                pf.map(pf, setSizeFieldMapper(pf.patch().size()));
            }
        }
    }
}


template<template<class> class GeoField>
void Foam::fvMeshStitcher::resizePatchFields()
{
    #define ResizePatchFields(Type, nullArg) \
        resizePatchFields<Type, GeoField>();
    FOR_ALL_FIELD_TYPES(ResizePatchFields);
    #undef ResizePatchFields
}


template<class Type>
void Foam::fvMeshStitcher::preConformSurfaceFields()
{
    UPtrList<SurfaceField<Type>> fields(mesh_.curFields<SurfaceField<Type>>());

    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];

        for (label ti=0; ti<=field.nOldTimes(false); ti++)
        {
            conformedFvsPatchField<Type>::conform
            (
                boundaryFieldRefNoUpdate(field.oldTime(ti))
            );
        }
    }
}


template<class Type>
void Foam::fvMeshStitcher::preConformVolFields()
{
    UPtrList<VolField<Type>> fields(mesh_.curFields<VolField<Type>>());

    forAll(fields, i)
    {
        VolField<Type>& field = fields[i];

        for (label ti=0; ti<=field.nOldTimes(false); ti++)
        {
            conformedFvPatchField<Type>::conform
            (
                boundaryFieldRefNoUpdate(field.oldTime(ti))
            );
        }
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformSurfaceFields()
{
    UPtrList<SurfaceField<Type>> fields(mesh_.curFields<SurfaceField<Type>>());

    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];

        for (label ti=0; ti<=field.nOldTimes(false); ti++)
        {
            conformedFvsPatchField<Type>::unconform
            (
                boundaryFieldRefNoUpdate(field.oldTime(ti))
            );
        }
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformVolFields()
{
    UPtrList<VolField<Type>> fields(mesh_.curFields<VolField<Type>>());

    forAll(fields, i)
    {
        VolField<Type>& field = fields[i];

        for (label ti=0; ti<=field.nOldTimes(false); ti++)
        {
            conformedFvPatchField<Type>::unconform
            (
                boundaryFieldRefNoUpdate(field.oldTime(ti))
            );
        }
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformEvaluateVolFields()
{
    auto evaluate = [](const typename VolField<Type>::Patch& pf)
    {
        return
            (
                isA<nonConformalFvPatch>(pf.patch())
             && pf.type() == pf.patch().patch().type()
             && polyPatch::constraintType(pf.patch().patch().type())
            )
         || isA<nonConformalErrorFvPatch>(pf.patch());
    };

    UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        const label nReq = Pstream::nRequests();

        forAll(mesh_.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                boundaryFieldRefNoUpdate(fields[i])[patchi];

            if (evaluate(pf))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mesh_.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                boundaryFieldRefNoUpdate(fields[i])[patchi];

            if (evaluate(pf))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class GeoField>
typename GeoField::Boundary& Foam::fvMeshStitcher::boundaryFieldRefNoUpdate
(
    GeoField& fld
)
{
    return const_cast<typename GeoField::Boundary&>(fld.boundaryField());
}


// ************************************************************************* //
