/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianPatch.H"
#include "LagrangianPatchField.H"
#include "LagrangianMesh.H"
#include "LagrangianSubFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianPatchField<Type>::LagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    patch_(p),
    internalIo_(iIo),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


template<class Type>
Foam::LagrangianPatchField<Type>::LagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary&
)
:
    patch_(p),
    internalIo_(iIo),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


template<class Type>
Foam::LagrangianPatchField<Type>::LagrangianPatchField
(
    const LagrangianPatchField<Type>& ptf
)
:
    patch_(ptf.patch_),
    internalIo_(ptf.internalIo_),
    internalField_(ptf.internalField_),
    internalNonDynamicField_(ptf.internalNonDynamicField_)
{}


template<class Type>
Foam::LagrangianPatchField<Type>::LagrangianPatchField
(
    const LagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    patch_(ptf.patch_),
    internalIo_(iIo),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


template<class Type>
Foam::LagrangianPatchField<Type>::LagrangianPatchField
(
    const LagrangianPatchField<Type>& ptf,
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    patch_(p),
    internalIo_(iIo),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::LagrangianPatchField<Type>>
Foam::LagrangianPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const LagrangianPatch& p,
    const regIOobject& iIo
)
{
    if (debug)
    {
        InfoInFunction << "Constructing LagrangianPatchField<Type>" << endl;
    }

    typename LagrangianPatchConstructorTable::iterator cstrIter =
        LagrangianPatchConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == LagrangianPatchConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown patchField type "
            << patchFieldType << nl << nl
            << "Valid patchField types are :" << endl
            << LagrangianPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        typename LagrangianPatchConstructorTable::iterator patchTypeCstrIter =
            LagrangianPatchConstructorTablePtr_->find(p.type());

        if (patchTypeCstrIter != LagrangianPatchConstructorTablePtr_->end())
        {
            return patchTypeCstrIter()(p, iIo);
        }
        else
        {
            return cstrIter()(p, iIo);
        }
    }
    else
    {
        return cstrIter()(p, iIo);
    }
}


template<class Type>
Foam::autoPtr<Foam::LagrangianPatchField<Type>>
Foam::LagrangianPatchField<Type>::New
(
    const word& patchFieldType,
    const LagrangianPatch& p,
    const regIOobject& iIo
)
{
    return New(patchFieldType, word::null, p, iIo);
}


template<class Type>
Foam::autoPtr<Foam::LagrangianPatchField<Type>>
Foam::LagrangianPatchField<Type>::New
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvsPatchField<Type>" << endl;
    }

    const word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowGenericLagrangianPatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("generic");
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    if
    (
        !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        typename dictionaryConstructorTable::iterator patchTypeCstrIter
            = dictionaryConstructorTablePtr_->find(p.type());

        if
        (
            patchTypeCstrIter != dictionaryConstructorTablePtr_->end()
         && patchTypeCstrIter() != cstrIter()
        )
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "inconsistent patch and patchField types for \n"
                   "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    return cstrIter()(p, iIo, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianPatchField<Type>::~LagrangianPatchField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::LagrangianPatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh();
}


template<class Type>
const Foam::LagrangianPatch& Foam::LagrangianPatchField<Type>::patch() const
{
    return patch_;
}


template<class Type>
const Foam::dimensionSet&
Foam::LagrangianPatchField<Type>::internalDimensions() const
{
    if (notNull(internalField_))
    {
        return internalField_.dimensions();
    }

    if (notNull(internalNonDynamicField_))
    {
        return internalNonDynamicField_.dimensions();
    }

    FatalErrorInFunction
        << "Dimensions of internal object " << internalIo_.name()
        << " could not be determined"
        << exit(FatalError);

    return NullObjectRef<dimensionSet>();
}


template<class Type>
const Foam::LagrangianInternalDynamicField<Type>&
Foam::LagrangianPatchField<Type>::internalField() const
{
    if (notNull(internalField_))
    {
        return internalField_;
    }

    FatalErrorInFunction
        << "Internal field " << internalIo_.name() << " is not of type "
        << LagrangianInternalDynamicField<Type>::typeName
        << exit(FatalError);

    return NullObjectRef<LagrangianInternalDynamicField<Type>>();
}


template<class Type>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::LagrangianPatchField<Type>::internalSubField() const
{
    if (notNull(internalField_))
    {
        return patch_.mesh().sub(internalField_);
    }

    if (notNull(internalNonDynamicField_))
    {
        return patch_.mesh().sub(internalNonDynamicField_);
    }

    FatalErrorInFunction
        << "Internal field " << internalIo_.name() << " is not of type "
        << LagrangianInternalDynamicField<Type>::typeName
        << exit(FatalError);

    return tmp<LagrangianSubSubField<Type>>(nullptr);
}


template<class Type>
Foam::SubField<Type> Foam::LagrangianPatchField<Type>::primitiveSubField() const
{
    if (notNull(internalField_))
    {
        return patch_.mesh().sub(internalField_.primitiveField());
    }

    if (notNull(internalNonDynamicField_))
    {
        return patch_.mesh().sub(internalNonDynamicField_.primitiveField());
    }

    FatalErrorInFunction
        << "Internal field " << internalIo_.name() << " is not of type "
        << LagrangianInternalDynamicField<Type>::typeName
        << exit(FatalError);

    return SubField<Type>(NullObjectRef<UList<Type>>(), 0, 0);
}


template<class Type>
void Foam::LagrangianPatchField<Type>::reset(const LagrangianPatchField<Type>&)
{}


template<class Type>
void Foam::LagrangianPatchField<Type>::initEvaluate
(
    PstreamBuffers&,
    const LagrangianScalarInternalDynamicField& fraction
)
{}


template<class Type>
void Foam::LagrangianPatchField<Type>::evaluate
(
    PstreamBuffers&,
    const LagrangianScalarInternalDynamicField& fraction
)
{}


template<class Type>
void Foam::LagrangianPatchField<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", type());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::LagrangianPatchField<Type>::operator()() const
{
    return internalSubField();
}


#define MEMBER_OPERATOR(MemberOp, FieldOp, OtherType, OtherDimensions)         \
                                                                               \
template<class Type>                                                           \
void Foam::LagrangianPatchField<Type>::operator MemberOp                       \
(                                                                              \
    const LagrangianPatchField<OtherType>& ptf                                 \
)                                                                              \
{                                                                              \
    if (&patch_ != &(ptf.patch_))                                              \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "different patches for LagrangianPatchField<Type>s"             \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    if (!patch_.boundaryMesh().mesh().changing()) return;                      \
                                                                               \
    OtherDimensions == ptf.internalDimensions();                               \
                                                                               \
    LagrangianSubSubField<Type> thisSsf(internalSubField());                   \
    thisSsf FieldOp ptf.internalSubField();                                    \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::LagrangianPatchField<Type>::operator MemberOp                       \
(                                                                              \
    const LagrangianSubField<OtherType>& sf                                    \
)                                                                              \
{                                                                              \
    if (!patch_.boundaryMesh().mesh().changing()) return;                      \
                                                                               \
    OtherDimensions == sf.dimensions();                                        \
                                                                               \
    SubField<Type> thisSf(primitiveSubField());                                \
    thisSf FieldOp sf.primitiveField();                                        \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::LagrangianPatchField<Type>::operator MemberOp                       \
(                                                                              \
    const LagrangianSubSubField<OtherType>& ssf                                \
)                                                                              \
{                                                                              \
    if (!patch_.boundaryMesh().mesh().changing()) return;                      \
                                                                               \
    OtherDimensions == ssf.dimensions();                                       \
                                                                               \
    SubField<Type> thisSf(primitiveSubField());                                \
    thisSf FieldOp ssf.primitiveField();                                       \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::LagrangianPatchField<Type>::operator MemberOp                       \
(                                                                              \
    const UList<OtherType>& f                                                  \
)                                                                              \
{                                                                              \
    if (!patch_.boundaryMesh().mesh().changing()) return;                      \
                                                                               \
    SubField<Type> thisSf(primitiveSubField());                                \
    thisSf FieldOp f;                                                          \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::LagrangianPatchField<Type>::operator MemberOp(const OtherType& t)   \
{                                                                              \
    if (!patch_.boundaryMesh().mesh().changing()) return;                      \
                                                                               \
    SubField<Type> thisSf(primitiveSubField());                                \
    thisSf FieldOp t;                                                          \
}

MEMBER_OPERATOR(=, =, Type, internalDimensions());
MEMBER_OPERATOR(==, =, Type, internalDimensions());
MEMBER_OPERATOR(+=, +=, Type, internalDimensions());
MEMBER_OPERATOR(-=, -=, Type, internalDimensions());
MEMBER_OPERATOR(*=, *=, scalar, dimless);
MEMBER_OPERATOR(/=, /=, scalar, dimless);

#undef MEMBER_OPERATOR


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LagrangianPatchField<Type>& ptf
)
{
    ptf.write(os);

    os.check
    (
        "Ostream& operator<<(Ostream&, const LagrangianPatchField<Type>&)"
    );

    return os;
}


// ************************************************************************* //
