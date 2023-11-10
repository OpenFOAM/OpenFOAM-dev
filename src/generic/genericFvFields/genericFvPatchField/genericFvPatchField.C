/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "genericFvPatchField.H"
#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bool insertUniformTypeField
(
    const scalarList& components,
    const word& keyword,
    const label size,
    HashPtrTable<Field<Type>>& typeFields
)
{
    if (components.size() != Type::nComponents)
    {
        return false;
    }

    Type t;
    forAll(t, i)
    {
        t[i] = components[i];
    }

    typeFields.insert(keyword, new Field<Type>(size, t));

    return true;
}


template<>
bool insertUniformTypeField
(
    const scalarList& components,
    const word& keyword,
    const label size,
    HashPtrTable<Field<scalar>>& typeFields
)
{
    return false;
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    genericFieldBase(dict.lookup("type")),
    calculatedFvPatchField<Type>(p, iF, dict),
    dict_(dict)
{
    if (!dict.found("value"))
    {
        FatalIOErrorInFunction(dict)
            << "\n    Cannot find 'value' entry"
            << " on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl
            << "    which is required to set the"
               " values of the generic patch field." << nl
            << "    (Actual type " << actualTypeName() << ")" << nl
            << "\n    Please add the 'value' entry to the write function "
               "of the user-defined boundary-condition\n"
            << exit(FatalIOError);
    }

    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type" && iter().keyword() != "value")
        {
            if
            (
                iter().isStream()
             && iter().stream().size()
            )
            {
                ITstream& is = iter().stream();

                // Read first token
                token firstToken(is);

                if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "nonuniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isCompound())
                    {
                        if
                        (
                            fieldToken.isLabel()
                         && fieldToken.labelToken() == 0
                        )
                        {
                            scalarFields_.insert
                            (
                                iter().keyword(),
                                new scalarField(0)
                            );
                        }
                        else
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    token following 'nonuniform' "
                                  "is not a compound"
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                            << exit(FatalIOError);
                        }
                    }

                    #define ReadTypeField(Type, nullArg)                       \
                        else if                                                \
                        (                                                      \
                            fieldToken.compoundToken().type()                  \
                         == token::Compound<List<Type>>::typeName              \
                        )                                                      \
                        {                                                      \
                            Field<Type>* fPtr = new Field<Type>;               \
                            fPtr->transfer                                     \
                            (                                                  \
                                dynamicCast<token::Compound<List<Type>>>       \
                                (                                              \
                                    fieldToken.transferCompoundToken(is)       \
                                )                                              \
                            );                                                 \
                                                                               \
                            if (fPtr->size() != this->size())                  \
                            {                                                  \
                                FatalIOErrorInFunction(dict)                   \
                                    << "\n    size of field "                  \
                                    << iter().keyword()                        \
                                    << " (" << fPtr->size() << ')'             \
                                    << " is not the same size as the patch ("  \
                                    << this->size() << ')'                     \
                                    << "\n    on patch "                       \
                                    << this->patch().name()                    \
                                    << " of field "                            \
                                    << this->internalField().name()            \
                                    << " in file "                             \
                                    << this->internalField().objectPath()      \
                                    << exit(FatalIOError);                     \
                            }                                                  \
                                                                               \
                            Type##Fields_.insert(iter().keyword(), fPtr);      \
                        }
                    FOR_ALL_FIELD_TYPES(ReadTypeField)
                    #undef ReadTypeField

                    else
                    {
                        FatalIOErrorInFunction(dict)
                            << "\n    compound " << fieldToken.compoundToken()
                            << " not supported"
                            << "\n    on patch " << this->patch().name()
                            << " of field "
                            << this->internalField().name()
                            << " in file "
                            << this->internalField().objectPath()
                            << exit(FatalIOError);
                    }
                }
                else if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "uniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isPunctuation())
                    {
                        scalarFields_.insert
                        (
                            iter().keyword(),
                            new scalarField
                            (
                                this->size(),
                                fieldToken.number()
                            )
                        );
                    }
                    else
                    {
                        // Read as a scalarList
                        is.putBack(fieldToken);
                        scalarList l(is);

                        #define InsertUniformTypeField(Type, nullArg)          \
                         || insertUniformTypeField                             \
                            (                                                  \
                                l,                                             \
                                iter().keyword(),                              \
                                this->size(),                                  \
                                Type##Fields_                                  \
                            )
                        if (!(0 FOR_ALL_FIELD_TYPES(InsertUniformTypeField)))
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    unrecognised native type " << l
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }
                        #undef InsertUniformTypeField
                    }
                }
            }
        }
    }
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const genericFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    genericFieldBase(ptf),
    calculatedFvPatchField<Type>(ptf, p, iF, mapper),
    dict_(ptf.dict_)
{
    #define MapTypeFields(Type, nullArg)                                       \
        forAllConstIter(HashPtrTable<Field<Type>>, ptf.Type##Fields_, iter)    \
        {                                                                      \
            Type##Fields_.insert                                               \
            (                                                                  \
                iter.key(),                                                    \
                mapper(*iter()).ptr()                                          \
            );                                                                 \
        }
    FOR_ALL_FIELD_TYPES(MapTypeFields);
    #undef MapTypeFields
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const genericFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    genericFieldBase(ptf),
    calculatedFvPatchField<Type>(ptf, iF),
    dict_(ptf.dict_),
    scalarFields_(ptf.scalarFields_),
    vectorFields_(ptf.vectorFields_),
    sphericalTensorFields_(ptf.sphericalTensorFields_),
    symmTensorFields_(ptf.symmTensorFields_),
    tensorFields_(ptf.tensorFields_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    calculatedFvPatchField<Type>::map(ptf, mapper);

    const genericFvPatchField<Type>& dptf =
        refCast<const genericFvPatchField<Type>>(ptf);

    #define MapTypeFields(Type, nullArg)                                       \
        forAllIter(HashPtrTable<Field<Type>>, Type##Fields_, iter)             \
        {                                                                      \
            HashPtrTable<Field<Type>>::const_iterator dptfIter =               \
                dptf.Type##Fields_.find(iter.key());                           \
                                                                               \
            if (dptfIter != dptf.Type##Fields_.end())                          \
            {                                                                  \
                mapper(*iter(), *dptfIter());                                  \
            }                                                                  \
        }
    FOR_ALL_FIELD_TYPES(MapTypeFields);
    #undef MapTypeFields
}


template<class Type>
void Foam::genericFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    calculatedFvPatchField<Type>::reset(ptf);

    const genericFvPatchField<Type>& dptf =
        refCast<const genericFvPatchField<Type>>(ptf);

    #define ResetTypeFields(Type, nullArg)                                     \
        forAllIter(HashPtrTable<Field<Type>>, Type##Fields_, iter)             \
        {                                                                      \
            HashPtrTable<Field<Type>>::const_iterator dptfIter =               \
                dptf.Type##Fields_.find(iter.key());                           \
                                                                               \
            if (dptfIter != dptf.Type##Fields_.end())                          \
            {                                                                  \
                iter()->reset(*dptfIter());                                    \
            }                                                                  \
        }
    FOR_ALL_FIELD_TYPES(ResetTypeFields);
    #undef ResetTypeFields
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName() << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName() << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName() << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName() << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
void Foam::genericFvPatchField<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", actualTypeName());

    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type" && iter().keyword() != "value")
        {
            if
            (
                iter().isStream()
             && iter().stream().size()
             && iter().stream()[0].isWord()
             && iter().stream()[0].wordToken() == "nonuniform"
            )
            {
                #define WriteTypeFieldEntry(Type, nullArg)                     \
                    else if (Type##Fields_.found(iter().keyword()))            \
                    {                                                          \
                        writeEntry                                             \
                        (                                                      \
                            os,                                                \
                            iter().keyword(),                                  \
                            *Type##Fields_.find(iter().keyword())()            \
                        );                                                     \
                    }
                if (false) {} FOR_ALL_FIELD_TYPES(WriteTypeFieldEntry)
                #undef WriteTypeFieldEntry
            }
            else
            {
               iter().write(os);
            }
        }
    }

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
