/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "GeometricField.H"
#include "Time.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "localIOdictionary.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define checkFieldAssignment(gf1, gf2)                                         \
                                                                               \
    if                                                                         \
    (                                                                          \
        static_cast<const regIOobject*>(&gf1)                                  \
     == static_cast<const regIOobject*>(&gf2)                                  \
    )                                                                          \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "attempted assignment to self for field "                       \
            << (gf1).name() << abort(FatalError);                              \
    }


#define checkFieldOperation(gf1, gf2, op)                                      \
                                                                               \
    if ((gf1).mesh() != (gf2).mesh())                                          \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "different mesh for fields "                                    \
            << (gf1).name() << " and " << (gf2).name()                         \
            << " during operation " <<  op                                     \
            << abort(FatalError);                                              \
    }


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::readFields
(
    const dictionary& dict
)
{
    Internal::readField(dict, "internalField");

    boundaryField_.readField(*this, dict.subDict("boundaryField"));

    // Don't use subOrEmptyDict here, or line numbers will be lost from any IO
    // error messages. Use an actual sub-dict reference here if possible.
    if (dict.found("sources"))
    {
        sources_.readField(*this, dict.subDict("sources"));
    }
    else
    {
        sources_.readField(*this, dictionary(dict.name()/"sources", dict));
    }

    if (dict.found("referenceLevel"))
    {
        Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));

        Field<Type>::operator+=(fieldAverage);

        forAll(boundaryField_, patchi)
        {
            boundaryField_[patchi] == boundaryField_[patchi] + fieldAverage;
        }
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::readFields()
{
    const localIOdictionary dict
    (
        IOobject
        (
            this->name(),
            this->instance(),
            this->local(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        typeName
    );

    this->close();

    readFields(dict);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
bool Foam::GeometricField<Type, GeoMesh, PrimitiveField>::readIfPresent()
{
    if
    (
        this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        WarningInFunction
            << "read option IOobject::MUST_READ or MUST_READ_IF_MODIFIED"
            << " suggests that a read constructor for field " << this->name()
            << " would be more appropriate." << endl;
    }
    else if
    (
        this->readOpt() == IOobject::READ_IF_PRESENT
     && this->headerOk()
    )
    {
        readFields();

        // Check compatibility between field and mesh
        if (this->size() != GeoMesh::size(this->mesh()))
        {
            FatalIOErrorInFunction(this->readStream(typeName))
                << "   number of field elements = " << this->size()
                << " number of mesh elements = "
                << GeoMesh::size(this->mesh())
                << exit(FatalIOError);
        }

        readOldTimeIfPresent();

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType
)
:
    Internal(io, mesh, ds, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType),
    sources_()
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
:
    Internal(io, mesh, ds, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes),
    sources_(*this, fieldSourceTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
:
    Internal(io, mesh, dt, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType),
    sources_()
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
:
    Internal(io, mesh, dt, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes),
    sources_(*this, fieldSourceTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Internal& diField,
    const PtrList<Patch>& ptfl,
    const HashPtrTable<Source>& stft
)
:
    Internal(io, diField, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, ptfl),
    sources_(*this, stft)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const PrimitiveField<Type>& iField,
    const PtrList<Patch>& ptfl,
    const HashPtrTable<Source>& stft
)
:
    Internal(io, mesh, ds, iField),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, ptfl),
    sources_(*this, stft)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh
)
:
    Internal(io, mesh, dimless, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary()),
    sources_()
{
    readFields();

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorInFunction(this->readStream(typeName))
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        InfoInFunction
            << "Finishing read-construction of" << endl << this->info() << endl;
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& dict
)
:
    Internal(io, mesh, dimless, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary()),
    sources_()
{
    readFields(dict);

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalErrorInFunction
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    if (debug)
    {
        InfoInFunction
            << "Finishing dictionary-construct of "
            << endl << this->info() << endl;
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
)
:
    Internal(gf),
    OldTimeField<GeometricField>(gf),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_),
    sources_(*this, gf.sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy" << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
:
    Internal(gf),
    OldTimeField<GeometricField>(gf),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_),
    sources_(*this, gf.sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy" << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    GeometricField<Type, GeoMesh, PrimitiveField>&& gf
)
:
    Internal(move(gf)),
    OldTimeField<GeometricField>(move(gf)),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_),
    sources_(*this, gf.sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing by moving" << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
:
    Internal
    (
        const_cast<GeometricField<Type, GeoMesh, PrimitiveField>&>(tgf()),
        tgf.isTmp()
    ),
    OldTimeField<GeometricField>(tgf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_),
    sources_(*this, tgf().sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp" << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
:
    Internal(io, gf, false),
    OldTimeField<GeometricField>(gf.timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_),
    sources_(*this, gf.sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    if (!readIfPresent())
    {
        copyOldTimes(io, gf);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, GeoMesh, PrimitiveField>&>(tgf()),
        tgf.isTmp(),
        false
    ),
    OldTimeField<GeometricField>(tgf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_),
    sources_(*this, tgf().sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting IO params"
            << endl << this->info() << endl;
    }

    tgf.clear();

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const word& newName,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
:
    Internal(newName, gf),
    OldTimeField<GeometricField>(gf.timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_),
    sources_(*this, gf.sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting name"
            << endl << this->info() << endl;
    }

    copyOldTimes(newName, gf);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const word& newName,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
:
    Internal
    (
        newName,
        const_cast<GeometricField<Type, GeoMesh, PrimitiveField>&>(tgf()),
        tgf.isTmp()
    ),
    OldTimeField<GeometricField>(tgf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_),
    sources_(*this, tgf().sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting name"
            << endl << this->info() << endl;
    }

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf,
    const word& patchFieldType
)
:
    Internal(io, gf, false),
    OldTimeField<GeometricField>(gf.timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType),
    sources_(*this, gf.sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent())
    {
        copyOldTimes(io, gf);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf,
    const word& patchFieldType
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, GeoMesh, PrimitiveField>&>(tgf()),
        tgf.isTmp(),
        false
    ),
    OldTimeField<GeometricField>(tgf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType),
    sources_(*this, tgf().sources_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    boundaryField_ == tgf().boundaryField_;

    tgf.clear();

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const word& patchFieldType
)
:
    Internal(io, df, false),
    OldTimeField<GeometricField>(this->time().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType),
    sources_()
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    if (!readIfPresent())
    {
        boundaryField_.evaluate();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const tmp<Internal>& tdf,
    const word& patchFieldType
)
:
    Internal(io, const_cast<Internal&>(tdf()), tdf.isTmp(), false),
    OldTimeField<GeometricField>(tdf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType),
    sources_()
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    if (!readIfPresent())
    {
        boundaryField_.evaluate();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
:
    Internal(io, gf, false),
    OldTimeField<GeometricField>(gf.timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    ),
    sources_(*this, fieldSourceTypes)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params and patch types"
            << endl << this->info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent())
    {
        copyOldTimes(io, gf);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, GeoMesh, PrimitiveField>&>(tgf()),
        tgf.isTmp(),
        false
    ),
    OldTimeField<GeometricField>(tgf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    ),
    sources_(*this, fieldSourceTypes)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting IO params and patch types"
            << endl << this->info() << endl;
    }

    boundaryField_ == tgf().boundaryField_;

    tgf.clear();

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
:
    Internal(io, df, false),
    OldTimeField<GeometricField>(df.timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    ),
    sources_(*this, fieldSourceTypes)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from internal field and patch types"
            << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::GeometricField
(
    const IOobject& io,
    const tmp<Internal>& tdf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
:
    Internal(io, const_cast<Internal&>(tdf()), tdf.isTmp(), false),
    OldTimeField<GeometricField>(tdf().timeIndex()),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    ),
    sources_(*this, fieldSourceTypes)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp internal field and patch types"
            << endl << this->info() << endl;
    }

    tdf.clear();

    readIfPresent();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::clone() const
{
    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>(*this)
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::cloneUnSliced() const
{
    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                this->name(),
                this->mesh().thisDb().time().name(),
                this->mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            Patch::calculatedType()
        )
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Internal& diField,
    const PtrList<Patch>& ptfl,
    const HashPtrTable<Source>& stft
)
{
    const bool cacheTmp = diField.mesh().thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                diField.mesh().thisDb().time().name(),
                diField.mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            diField,
            ptfl,
            stft
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().name(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            mesh,
            ds,
            patchFieldType
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().name(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            mesh,
            dt,
            patchFieldType
        ),
        cacheTmp
    );
}



template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().name(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            mesh,
            dt,
            patchFieldTypes,
            actualPatchTypes,
            fieldSourceTypes
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    const bool cacheTmp = tgf().db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                tgf().instance(),
                tgf().local(),
                tgf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            tgf
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf,
    const word& patchFieldType
)
{
    const bool cacheTmp = gf.db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                gf.instance(),
                gf.local(),
                gf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            gf,
            patchFieldType
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf,
    const word& patchFieldType
)
{
    const bool cacheTmp = tgf().db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                tgf().instance(),
                tgf().local(),
                tgf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            tgf,
            patchFieldType
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const word& patchFieldType
)
{
    const bool cacheTmp = df.db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                df.instance(),
                df.local(),
                df.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            df,
            patchFieldType
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const tmp<Internal>& tdf,
    const word& patchFieldType
)
{
    const bool cacheTmp = tdf().db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                tdf().instance(),
                tdf().local(),
                tdf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            tdf,
            patchFieldType
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
{
    const bool cacheTmp = gf.db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                gf.instance(),
                gf.local(),
                gf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            gf,
            patchFieldTypes,
            actualPatchTypes,
            fieldSourceTypes
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
{
    const bool cacheTmp = tgf().db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                tgf().instance(),
                tgf().local(),
                tgf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            tgf,
            patchFieldTypes,
            actualPatchTypes,
            fieldSourceTypes
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
{
    const bool cacheTmp = df.db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                df.instance(),
                df.local(),
                df.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            df,
            patchFieldTypes,
            actualPatchTypes,
            fieldSourceTypes
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const tmp<Internal>& tdf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes,
    const HashTable<word>& fieldSourceTypes
)
{
    const bool cacheTmp = tdf().db().cacheTemporaryObject(newName);

    return tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
    (
        new GeometricField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                tdf().instance(),
                tdf().local(),
                tdf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            tdf,
            patchFieldTypes,
            actualPatchTypes,
            fieldSourceTypes
        ),
        cacheTmp
    );
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::~GeometricField()
{
    this->db().cacheTemporaryObject(*this);

    clearPrevIter();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
typename Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Internal&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::internalFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Internal::FieldType&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::primitiveFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Boundary&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::boundaryFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return boundaryField_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Boundary&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::
boundaryFieldRefNoStoreOldTimes()
{
    this->setUpToDate();
    return boundaryField_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Sources&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::sourcesRef()
{
    this->setUpToDate();
    storeOldTimes();
    return sources_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::storePrevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        if (debug)
        {
            InfoInFunction
                << "Allocating previous iteration field" << endl
                << this->info() << endl;
        }

        fieldPrevIterPtr_ =
            new GeometricField<Type, GeoMesh, Field>
            (
                this->name() + "PrevIter",
                *this
            );
    }
    else
    {
        *fieldPrevIterPtr_ == *this;
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::clearPrevIter()
{
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
const Foam::GeometricField<Type, GeoMesh, Foam::Field>&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::prevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorInFunction
            << "previous iteration field" << endl << this->info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::
correctBoundaryConditions()
{
    this->setUpToDate();
    storeOldTimes();
    boundaryField_.evaluate();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::reset
(
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
{
    Internal::reset(gf);

    boundaryField_.reset(gf.boundaryField());
    sources_.reset(*this, gf.sources());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::reset
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    checkFieldAssignment(*this, gf);

    this->dimensions() = gf.dimensions();

    if (tgf.isTmp())
    {
        PrimitiveField<Type>::transfer(tgf.ref());
    }
    else
    {
        PrimitiveField<Type>::operator=(gf.primitiveField());
    }

    boundaryField_.reset(gf.boundaryField());
    sources_.reset(*this, gf.sources());

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::reset
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>& tgf
)
{
    reset(tgf());

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
bool Foam::GeometricField<Type, GeoMesh, PrimitiveField>::needReference() const
{
    // Search all boundary conditions, if any are
    // fixed-value or mixed (Robin) do not set reference level for solution.

    bool needRef = true;

    forAll(boundaryField_, patchi)
    {
        if (boundaryField_[patchi].fixesValue())
        {
            needRef = false;
            break;
        }
    }

    reduce(needRef, andOp<bool>());

    return needRef;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::relax
(
    const scalar alpha
)
{
    if (alpha < 1)
    {
        if (debug)
        {
            InfoInFunction
                << "Relaxing" << endl << this->info()
                << " by " << alpha << endl;
        }

        operator==(prevIter() + alpha*(*this - prevIter()));
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::scalar
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::relaxationFactor() const
{
    if
    (
        solutionControl::finalIteration(this->mesh())
     && this->mesh().solution().relaxField(this->name() + "Final")
    )
    {
        return this->mesh().solution().fieldRelaxationFactor
        (
            this->name() + "Final"
        );
    }
    else if (this->mesh().solution().relaxField(this->name()))
    {
        return this->mesh().solution().fieldRelaxationFactor(this->name());
    }
    else
    {
        return 1;
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::relax()
{
    relax(relaxationFactor());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::relax
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>& tgf,
    const scalar alpha
)
{
    if (alpha < 1)
    {
        if (debug)
        {
            InfoInFunction
                << "Relaxing" << endl << this->info()
                << " by " << alpha << endl;
        }

        operator==(*this + alpha*(tgf - *this));
    }
    else
    {
        operator==(tgf);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::relax
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>& tgf
)
{
    relax(tgf, relaxationFactor());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::word Foam::GeometricField<Type, GeoMesh, PrimitiveField>::select
(
    bool final
) const
{
    if (final)
    {
        return this->name() + "Final";
    }
    else
    {
        return this->name();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::writeMinMax
(
    Ostream& os
) const
{
    os  << "min/max(" << this->name() << ") = "
        << Foam::min(*this).value() << ", "
        << Foam::max(*this).value()
        << endl;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
bool Foam::GeometricField<Type, GeoMesh, PrimitiveField>::writeData
(
    Ostream& os
) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh, Foam::Field>>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::T() const
{
    tmp<GeometricField<Type, GeoMesh, Field>> result
    (
        GeometricField<Type, GeoMesh, Field>::New
        (
            this->name() + ".T()",
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::T(result.ref().primitiveFieldRef(), primitiveField());
    Foam::T(result.ref().boundaryFieldRef(), boundaryField());

    return result;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::GeometricField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        Foam::Field
    >
>
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::component
(
    const direction d
) const
{
    tmp<GeometricField<cmptType, GeoMesh, Field>> Component
    (
        GeometricField<cmptType, GeoMesh, Field>::New
        (
            this->name() + ".component(" + Foam::name(d) + ')',
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::component(Component.ref().primitiveFieldRef(), primitiveField(), d);
    Foam::component(Component.ref().boundaryFieldRef(), boundaryField(), d);

    return Component;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::replace
(
    const direction d,
    const GeometricField
    <
        typename GeometricField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        PrimitiveField2
     >& gcf
)
{
    primitiveFieldRef().replace(d, gcf.primitiveField());
    boundaryFieldRef().replace(d, gcf.boundaryField());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::replace
(
    const direction d,
    const tmp
    <
        GeometricField
        <
            typename GeometricField<Type, GeoMesh, PrimitiveField>::cmptType,
            GeoMesh,
            PrimitiveField2
        >
    >& gcf
)
{
    replace(d, gcf());
    gcf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::replace
(
    const direction d,
    const dimensioned<cmptType>& ds
)
{
    primitiveFieldRef().replace(d, ds.value());
    boundaryFieldRef().replace(d, ds.value());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::max
(
    const dimensioned<Type>& dt
)
{
    Foam::max(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::max(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::min
(
    const dimensioned<Type>& dt
)
{
    Foam::min(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::min(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::maxMin
(
    const dimensioned<Type>& minDt,
    const dimensioned<Type>& maxDt
)
{
    Foam::max(primitiveFieldRef(), primitiveField(), minDt.value());
    Foam::max(boundaryFieldRef(), boundaryField(), minDt.value());
    Foam::min(primitiveFieldRef(), primitiveField(), maxDt.value());
    Foam::min(boundaryFieldRef(), boundaryField(), maxDt.value());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::negate()
{
    primitiveFieldRef().negate();
    boundaryFieldRef().negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
)
{
    checkFieldAssignment(*this, gf);
    checkFieldOperation(*this, gf, "=");

    internalFieldRef() = gf.internalField();
    boundaryFieldRef() = gf.boundaryField();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    GeometricField<Type, GeoMesh, PrimitiveField>&& gf
)
{
    checkFieldAssignment(*this, gf);
    checkFieldOperation(*this, gf, "=");

    internalFieldRef() = move(gf.internalField());
    boundaryFieldRef() = move(gf.boundaryField());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
{
    checkFieldOperation(*this, gf, "=");

    internalFieldRef() = gf.internalField();
    boundaryFieldRef() = gf.boundaryField();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    checkFieldAssignment(*this, gf);
    checkFieldOperation(*this, gf, "=");

    this->dimensions() = gf.dimensions();

    if (tgf.isTmp())
    {
        primitiveFieldRef().transfer(tgf.ref());
    }
    else
    {
        primitiveFieldRef() = gf.primitiveField();
    }

    boundaryFieldRef() = gf.boundaryField();

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf = tgf();

    checkFieldOperation(*this, gf, "=");

    internalFieldRef() = gf.internalField();
    boundaryFieldRef() = gf.boundaryField();

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    const dimensioned<Type>& dt
)
{
    internalFieldRef() = dt;
    boundaryFieldRef() = dt.value();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator=
(
    const zero&
)
{
    internalFieldRef() = Zero;
    boundaryFieldRef() = Zero;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator==
(
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf
)
{
    checkFieldOperation(*this, gf, "==");

    internalFieldRef() = gf.internalField();
    boundaryFieldRef() == gf.boundaryField();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator==
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf = tgf();

    checkFieldOperation(*this, gf, "==");

    this->dimensions() = gf.dimensions();

    if (tgf.isTmp())
    {
        primitiveFieldRef().transfer(tgf.ref());
    }
    else
    {
        primitiveFieldRef() = gf.primitiveField();
    }

    boundaryFieldRef() == gf.boundaryField();

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator==
(
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>& tgf
)
{
    const GeometricField<Type, GeoMesh, PrimitiveField2>& gf = tgf();

    checkFieldOperation(*this, gf, "=");

    internalFieldRef() = gf.internalField();
    boundaryFieldRef() == gf.boundaryField();

    tgf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator==
(
    const dimensioned<Type>& dt
)
{
    internalFieldRef() = dt;
    boundaryFieldRef() == dt.value();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator==
(
    const zero&
)
{
    internalFieldRef() = Zero;
    boundaryFieldRef() == Zero;
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
template<template<class> class PrimitiveField2>                                \
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator op          \
(                                                                              \
    const GeometricField<TYPE, GeoMesh, PrimitiveField2>& gf                   \
)                                                                              \
{                                                                              \
    checkFieldOperation(*this, gf, #op);                                       \
                                                                               \
    internalFieldRef() op gf.internalField();                                  \
    boundaryFieldRef() op gf.boundaryField();                                  \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
template<template<class> class PrimitiveField2>                                \
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator op          \
(                                                                              \
    const tmp<GeometricField<TYPE, GeoMesh, PrimitiveField2>>& tgf             \
)                                                                              \
{                                                                              \
    operator op(tgf());                                                        \
    tgf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
void Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator op          \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    internalFieldRef() op dt;                                                  \
    boundaryFieldRef() op dt.value();                                          \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricField<Type, GeoMesh, PrimitiveField>& gf
)
{
    gf().writeData(os, "internalField");
    os  << nl;
    gf.boundaryField().writeEntry("boundaryField", os);

    if (!gf.sources_.empty())
    {
        os  << nl;
        gf.sources().writeEntry("sources", os);
    }

    // Check state of IOstream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const GeometricField<Type, GeoMesh, PrimitiveField>&)"
    );

    return (os);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf
)
{
    os << tgf();
    tgf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkFieldAssignment
#undef checkFieldOperation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricFieldFunctions.C"

// ************************************************************************* //
