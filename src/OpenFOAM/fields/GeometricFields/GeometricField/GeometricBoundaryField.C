/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "GeometricBoundaryField.H"
#include "GeometricFieldFwd.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "commSchedule.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::readField
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const dictionary& dict
)
{
    // Clear the boundary field if already initialised
    this->clear();

    this->setSize(bmesh_.size());

    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    // Construct a list of entry pointers for each patch
    UPtrList<const entry> patchEntries(this->size());

    // 1. Explicit patch names
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict() && !iter().keyword().isPattern())
        {
            const label patchi = bmesh_.findIndex(iter().keyword());

            if (patchi != -1)
            {
                patchEntries.set(patchi, &(iter()));
            }
        }
    }

    // 2. Patch-groups
    //    Note: This is done in reverse order of the entries in the dictionary,
    //    so that it is consistent with dictionary wildcard behaviour.
    if (dict.size())
    {
        for
        (
            IDLList<entry>::const_reverse_iterator iter = dict.rbegin();
            iter != dict.rend();
            ++iter
        )
        {
            if (iter().isDict() && !iter().keyword().isPattern())
            {
                const labelList patchIDs =
                    bmesh_.findIndices(wordRe(iter().keyword()), true);

                forAll(patchIDs, i)
                {
                    const label patchi = patchIDs[i];

                    if (!patchEntries.set(patchi))
                    {
                        patchEntries.set(patchi, &(iter()));
                    }
                }
            }
        }
    }

    // 3. Empty patches
    //    These take precedence over wildcards
    //    (... apparently. Why not wedges and/or other constraints too?)
    forAll(bmesh_, patchi)
    {
        if (!patchEntries.set(patchi))
        {
            if (bmesh_[patchi].type() == emptyPolyPatch::typeName)
            {
                patchEntries.set(patchi, NullObjectPtr<entry>());
            }
        }
    }

    // 4. Wildcards
    forAll(bmesh_, patchi)
    {
        if (!patchEntries.set(patchi))
        {
            const entry* ePtr =
                dict.lookupEntryPtr(bmesh_[patchi].name(), false, true);

            if (ePtr)
            {
                patchEntries.set(patchi, ePtr);
            }
        }
    }

    // Construct all the patches in order
    forAll(bmesh_, patchi)
    {
        if (patchEntries.set(patchi) && !isNull(patchEntries(patchi)))
        {
            this->set
            (
                patchi,
                Patch::New
                (
                    bmesh_[patchi],
                    field,
                    patchEntries[patchi].dict()
                )
            );
        }
        else if (patchEntries.set(patchi) && isNull(patchEntries[patchi]))
        {
            this->set
            (
                patchi,
                Patch::New
                (
                    emptyPolyPatch::typeName,
                    bmesh_[patchi],
                    field
                )
            );
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find patchField entry for "
                << bmesh_[patchi].name() << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh
)
:
    FieldField<GeoMesh::template PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const word& patchFieldType
)
:
    FieldField<GeoMesh::template PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    forAll(bmesh_, patchi)
    {
        this->set
        (
            patchi,
            Patch::New
            (
                patchFieldType,
                bmesh_[patchi],
                field
            )
        );
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const wordList& patchFieldTypes,
    const wordList& constraintTypes
)
:
    FieldField<GeoMesh::template PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    if
    (
        patchFieldTypes.size() != this->size()
     || (constraintTypes.size() && (constraintTypes.size() != this->size()))
    )
    {
        FatalErrorInFunction
            << "Incorrect number of patch type specifications given" << nl
            << "    Number of patches in mesh = " << bmesh.size()
            << " number of patch type specifications = "
            << patchFieldTypes.size()
            << abort(FatalError);
    }

    if (constraintTypes.size())
    {
        forAll(bmesh_, patchi)
        {
            this->set
            (
                patchi,
                Patch::New
                (
                    patchFieldTypes[patchi],
                    constraintTypes[patchi],
                    bmesh_[patchi],
                    field
                )
            );
        }
    }
    else
    {
        forAll(bmesh_, patchi)
        {
            this->set
            (
                patchi,
                Patch::New
                (
                    patchFieldTypes[patchi],
                    bmesh_[patchi],
                    field
                )
            );
        }
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const PtrList<Patch>& ptfl
)
:
    FieldField<GeoMesh::template PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    forAll(bmesh_, patchi)
    {
        this->set(patchi, ptfl[patchi].clone(field));
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const GeometricBoundaryField<Type, GeoMesh, PrimitiveField>& btf
)
:
    FieldField<GeoMesh::template PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    forAll(bmesh_, patchi)
    {
        this->set(patchi, btf[patchi].clone(field));
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const dictionary& dict
)
:
    FieldField<GeoMesh::template PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    readField(field, dict);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
GeometricBoundaryField
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field,
    const GeometricBoundaryField<Type, GeoMesh, PrimitiveField2>& btf
)
:
    FieldField<GeoMesh::template PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    forAll(bmesh_, patchi)
    {
        this->set(patchi, btf[patchi].clone(field));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::updateCoeffs()
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    forAll(*this, patchi)
    {
        this->operator[](patchi).updateCoeffs();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::evaluate()
{
    if (GeometricField<Type, GeoMesh, Field>::debug)
    {
        InfoInFunction << endl;
    }

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(*this, patchi)
        {
            this->operator[](patchi).initEvaluate(Pstream::defaultCommsType);
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

        forAll(*this, patchi)
        {
            this->operator[](patchi).evaluate(Pstream::defaultCommsType);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            if (patchSchedule[patchEvali].init)
            {
                this->operator[](patchSchedule[patchEvali].patch)
                    .initEvaluate(Pstream::commsTypes::scheduled);
            }
            else
            {
                this->operator[](patchSchedule[patchEvali].patch)
                    .evaluate(Pstream::commsTypes::scheduled);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::wordList
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::types() const
{
    const FieldField<GeoMesh::template PatchField, Type>& pff = *this;

    wordList Types(pff.size());

    forAll(pff, patchi)
    {
        Types[patchi] = pff[patchi].type();
    }

    return Types;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
boundaryInternalField() const
{
    typedef
        GeometricBoundaryField<Type, GeoMesh, PrimitiveField>
        resultType;

    tmp<resultType> tresult(new resultType(Internal::null(), *this));
    resultType& result = tresult.ref();

    forAll(*this, patchi)
    {
        result[patchi] == this->operator[](patchi).patchInternalField();
    }

    return tresult;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
boundaryNeighbourField() const
{
    typedef
        GeometricBoundaryField<Type, GeoMesh, PrimitiveField>
        resultType;

    tmp<resultType> tresult(new resultType(Internal::null(), *this));
    resultType& result = tresult.ref();

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        const label nReq = Pstream::nRequests();

        forAll(*this, patchi)
        {
            if (this->operator[](patchi).coupled())
            {
                this->operator[](patchi)
                    .initPatchNeighbourField(Pstream::defaultCommsType);
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

        forAll(*this, patchi)
        {
            if (this->operator[](patchi).coupled())
            {
                result[patchi] =
                    this->operator[](patchi)
                    .patchNeighbourField(Pstream::defaultCommsType);
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            if (this->operator[](patchSchedule[patchEvali].patch).coupled())
            {
                if (patchSchedule[patchEvali].init)
                {
                    this->operator[](patchSchedule[patchEvali].patch)
                        .initPatchNeighbourField(Pstream::defaultCommsType);
                }
                else
                {
                    result[patchSchedule[patchEvali].patch] =
                        this->operator[](patchSchedule[patchEvali].patch)
                        .patchNeighbourField(Pstream::defaultCommsType);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }

    return tresult;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::LduInterfaceFieldPtrsList<Type>
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::interfaces() const
{
    LduInterfaceFieldPtrsList<Type> interfaces(this->size());

    forAll(interfaces, patchi)
    {
        if (isA<LduInterfaceField<Type>>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
                &refCast<const LduInterfaceField<Type>>
                (
                    this->operator[](patchi)
                )
            );
        }
    }

    return interfaces;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::lduInterfaceFieldPtrsList
Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
scalarInterfaces() const
{
    lduInterfaceFieldPtrsList interfaces(this->size());

    forAll(interfaces, patchi)
    {
        if (isA<lduInterfaceField>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
                &refCast<const lduInterfaceField>
                (
                    this->operator[](patchi)
                )
            );
        }
    }

    return interfaces;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::reset
(
    const GeometricBoundaryField<Type, GeoMesh, PrimitiveField>& btf
)
{
    // Reset the number of patches in case the decomposition changed
    this->setSize(btf.size());

    const polyBoundaryMesh& pbm = bmesh_.mesh().mesh().boundaryMesh();

    forAll(*this, patchi)
    {
        // Construct new processor patch fields in case the decomposition
        // changed
        if (isA<processorPolyPatch>(pbm[patchi]))
        {
            this->set
            (
                patchi,
                btf[patchi].clone
                (
                    bmesh_[patchi],
                    this->operator[](0).internalField()
                )
            );
        }
        else
        {
            this->operator[](patchi).reset(btf[patchi]);
        }
    }
}



template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os  << keyword << nl << token::BEGIN_BLOCK << incrIndent << nl;

    forAll(*this, patchi)
    {
        os  << indent << this->operator[](patchi).patch().name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << this->operator[](patchi) << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_BLOCK << endl;

    // Check state of IOstream
    os.check
    (
        "GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::"
        "writeEntry(const word& keyword, Ostream& os) const"
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator=(const GeometricBoundaryField<Type, GeoMesh, PrimitiveField>& bf)
{
    FieldField<GeoMesh::template PatchField, Type>::operator=(bf);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator=(GeometricBoundaryField<Type, GeoMesh, PrimitiveField>&& bf)
{
    FieldField<GeoMesh::template PatchField, Type>::operator=(move(bf));
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator=(const FieldField<GeoMesh::template PatchField, Type>& ptff)
{
    FieldField<GeoMesh::template PatchField, Type>::operator=(ptff);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class OtherPatchField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator=(const FieldField<OtherPatchField, Type>& ptff)
{
    FieldField<GeoMesh::template PatchField, Type>::operator=(ptff);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator=(const Type& t)
{
    FieldField<GeoMesh::template PatchField, Type>::operator=(t);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator==(const GeometricBoundaryField<Type, GeoMesh, PrimitiveField>& bf)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator==(const FieldField<GeoMesh::template PatchField, Type>& ptff)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == ptff[patchi];
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class OtherPatchField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator==(const FieldField<OtherPatchField, Type>& ptff)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == ptff[patchi];
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void Foam::GeometricBoundaryField<Type, GeoMesh, PrimitiveField>::
operator==(const Type& t)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == t;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricBoundaryField<Type, GeoMesh, PrimitiveField>& bf
)
{
    os <<
        static_cast<const FieldField<GeoMesh::template PatchField, Type>&>(bf);
    return os;
}


// ************************************************************************* //
