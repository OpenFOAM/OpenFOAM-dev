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

#include "genericPointPatchField.H"
#include "fieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    genericFieldBase(dict.lookup("type")),
    calculatedPointPatchField<Type>(p, iF, dict),
    dict_(dict)
{
    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type")
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
                }
            }
        }
    }
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const fieldMapper& mapper
)
:
    genericFieldBase(ptf),
    calculatedPointPatchField<Type>(ptf, p, iF, mapper),
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
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    genericFieldBase(ptf),
    calculatedPointPatchField<Type>(ptf, iF),
    dict_(ptf.dict_)
    #define CopyTypeFields(Type, nullArg) \
        , Type##Fields_(ptf.Type##Fields_)
    FOR_ALL_FIELD_TYPES(CopyTypeFields)
    #undef CopyTypeFields
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericPointPatchField<Type>::map
(
    const pointPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    const genericPointPatchField<Type>& dptf =
        refCast<const genericPointPatchField<Type>>(ptf);

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
void Foam::genericPointPatchField<Type>::reset
(
    const pointPatchField<Type>& ptf
)
{
    const genericPointPatchField<Type>& dptf =
        refCast<const genericPointPatchField<Type>>(ptf);

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
void Foam::genericPointPatchField<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", actualTypeName());

    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type")
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
}


// ************************************************************************* //
