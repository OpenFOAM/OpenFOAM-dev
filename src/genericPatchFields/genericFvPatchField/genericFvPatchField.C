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

#include "genericFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(p, iF)
{
    FatalErrorInFunction
        << "Trying to construct an genericFvPatchField on patch "
        << this->patch().name()
        << " of field " << this->internalField().name()
        << abort(FatalError);
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    genericPatchField(dict.lookup("type")),
    calculatedFvPatchField<Type>(p, iF, dict),
    dict_(dict)
{
    if (!dict.found("value"))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "\n    Cannot find 'value' entry"
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
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    token following 'nonuniform' "
                                  "is not a compound"
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
                        fieldToken.compoundToken().type()
                     == token::Compound<List<scalar>>::typeName
                    )
                    {
                        scalarField* fPtr = new scalarField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<scalar>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        scalarFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<vector>>::typeName
                    )
                    {
                        vectorField* fPtr = new vectorField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<vector>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        vectorFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<sphericalTensor>>::typeName
                    )
                    {
                        sphericalTensorField* fPtr = new sphericalTensorField;
                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<sphericalTensor>>
                            >
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        sphericalTensorFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<symmTensor>>::typeName
                    )
                    {
                        symmTensorField* fPtr = new symmTensorField;
                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<symmTensor>>
                            >
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        symmTensorFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<tensor>>::typeName
                    )
                    {
                        tensorField* fPtr = new tensorField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<tensor>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        tensorFields_.insert(iter().keyword(), fPtr);
                    }
                    else
                    {
                        FatalIOErrorInFunction
                        (
                            dict
                        )   << "\n    compound " << fieldToken.compoundToken()
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
                        // Read as scalarList.
                        is.putBack(fieldToken);

                        scalarList l(is);

                        if (l.size() == vector::nComponents)
                        {
                            vector vs(l[0], l[1], l[2]);

                            vectorFields_.insert
                            (
                                iter().keyword(),
                                new vectorField(this->size(), vs)
                            );
                        }
                        else if (l.size() == sphericalTensor::nComponents)
                        {
                            sphericalTensor vs(l[0]);

                            sphericalTensorFields_.insert
                            (
                                iter().keyword(),
                                new sphericalTensorField(this->size(), vs)
                            );
                        }
                        else if (l.size() == symmTensor::nComponents)
                        {
                            symmTensor vs(l[0], l[1], l[2], l[3], l[4], l[5]);

                            symmTensorFields_.insert
                            (
                                iter().keyword(),
                                new symmTensorField(this->size(), vs)
                            );
                        }
                        else if (l.size() == tensor::nComponents)
                        {
                            tensor vs
                            (
                                l[0], l[1], l[2],
                                l[3], l[4], l[5],
                                l[6], l[7], l[8]
                            );

                            tensorFields_.insert
                            (
                                iter().keyword(),
                                new tensorField(this->size(), vs)
                            );
                        }
                        else
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    unrecognised native type " << l
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }
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
    const fvPatchFieldMapper& mapper
)
:
    genericPatchField(ptf),
    calculatedFvPatchField<Type>(ptf, p, iF, mapper),
    dict_(ptf.dict_)
{
    forAllConstIter
    (
        HashPtrTable<scalarField>,
        ptf.scalarFields_,
        iter
    )
    {
        scalarFields_.insert
        (
            iter.key(),
            mapper(*iter()).ptr()
        );
    }

    forAllConstIter
    (
        HashPtrTable<vectorField>,
        ptf.vectorFields_,
        iter
    )
    {
        vectorFields_.insert
        (
            iter.key(),
            mapper(*iter()).ptr()
        );
    }

    forAllConstIter
    (
        HashPtrTable<sphericalTensorField>,
        ptf.sphericalTensorFields_,
        iter
    )
    {
        sphericalTensorFields_.insert
        (
            iter.key(),
            mapper(*iter()).ptr()
        );
    }

    forAllConstIter
    (
        HashPtrTable<symmTensorField>,
        ptf.symmTensorFields_,
        iter
    )
    {
        symmTensorFields_.insert
        (
            iter.key(),
            mapper(*iter()).ptr()
        );
    }

    forAllConstIter
    (
        HashPtrTable<tensorField>,
        ptf.tensorFields_,
        iter
    )
    {
        tensorFields_.insert
        (
            iter.key(),
            mapper(*iter()).ptr()
        );
    }
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const genericFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    genericPatchField(ptf),
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
void Foam::genericFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    calculatedFvPatchField<Type>::autoMap(m);

    forAllIter
    (
        HashPtrTable<scalarField>,
        scalarFields_,
        iter
    )
    {
        m(*iter(), *iter());
    }

    forAllIter
    (
        HashPtrTable<vectorField>,
        vectorFields_,
        iter
    )
    {
        m(*iter(), *iter());
    }

    forAllIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields_,
        iter
    )
    {
        m(*iter(), *iter());
    }

    forAllIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields_,
        iter
    )
    {
        m(*iter(), *iter());
    }

    forAllIter
    (
        HashPtrTable<tensorField>,
        tensorFields_,
        iter
    )
    {
        m(*iter(), *iter());
    }
}


template<class Type>
void Foam::genericFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    calculatedFvPatchField<Type>::rmap(ptf, addr);

    const genericFvPatchField<Type>& dptf =
        refCast<const genericFvPatchField<Type>>(ptf);

    forAllIter
    (
        HashPtrTable<scalarField>,
        scalarFields_,
        iter
    )
    {
        HashPtrTable<scalarField>::const_iterator dptfIter =
            dptf.scalarFields_.find(iter.key());

        if (dptfIter != dptf.scalarFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<vectorField>,
        vectorFields_,
        iter
    )
    {
        HashPtrTable<vectorField>::const_iterator dptfIter =
            dptf.vectorFields_.find(iter.key());

        if (dptfIter != dptf.vectorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields_,
        iter
    )
    {
        HashPtrTable<sphericalTensorField>::const_iterator dptfIter =
            dptf.sphericalTensorFields_.find(iter.key());

        if (dptfIter != dptf.sphericalTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields_,
        iter
    )
    {
        HashPtrTable<symmTensorField>::const_iterator dptfIter =
            dptf.symmTensorFields_.find(iter.key());

        if (dptfIter != dptf.symmTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<tensorField>,
        tensorFields_,
        iter
    )
    {
        HashPtrTable<tensorField>::const_iterator dptfIter =
            dptf.tensorFields_.find(iter.key());

        if (dptfIter != dptf.tensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }
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

    forAllIter
    (
        HashPtrTable<scalarField>,
        scalarFields_,
        iter
    )
    {
        HashPtrTable<scalarField>::const_iterator dptfIter =
            dptf.scalarFields_.find(iter.key());

        if (dptfIter != dptf.scalarFields_.end())
        {
            iter()->reset(*dptfIter());
        }
    }

    forAllIter
    (
        HashPtrTable<vectorField>,
        vectorFields_,
        iter
    )
    {
        HashPtrTable<vectorField>::const_iterator dptfIter =
            dptf.vectorFields_.find(iter.key());

        if (dptfIter != dptf.vectorFields_.end())
        {
            iter()->reset(*dptfIter());
        }
    }

    forAllIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields_,
        iter
    )
    {
        HashPtrTable<sphericalTensorField>::const_iterator dptfIter =
            dptf.sphericalTensorFields_.find(iter.key());

        if (dptfIter != dptf.sphericalTensorFields_.end())
        {
            iter()->reset(*dptfIter());
        }
    }

    forAllIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields_,
        iter
    )
    {
        HashPtrTable<symmTensorField>::const_iterator dptfIter =
            dptf.symmTensorFields_.find(iter.key());

        if (dptfIter != dptf.symmTensorFields_.end())
        {
            iter()->reset(*dptfIter());
        }
    }

    forAllIter
    (
        HashPtrTable<tensorField>,
        tensorFields_,
        iter
    )
    {
        HashPtrTable<tensorField>::const_iterator dptfIter =
            dptf.tensorFields_.find(iter.key());

        if (dptfIter != dptf.tensorFields_.end())
        {
            iter()->reset(*dptfIter());
        }
    }
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
                if (scalarFields_.found(iter().keyword()))
                {
                    writeEntry
                    (
                        os,
                        iter().keyword(),
                        *scalarFields_.find(iter().keyword())()
                    );
                }
                else if (vectorFields_.found(iter().keyword()))
                {
                    writeEntry
                    (
                        os,
                        iter().keyword(),
                        *vectorFields_.find(iter().keyword())()
                    );
                }
                else if (sphericalTensorFields_.found(iter().keyword()))
                {
                    writeEntry
                    (
                        os,
                        iter().keyword(),
                        *sphericalTensorFields_.find(iter().keyword())()
                    );
                }
                else if (symmTensorFields_.found(iter().keyword()))
                {
                    writeEntry
                    (
                        os,
                        iter().keyword(),
                        *symmTensorFields_.find(iter().keyword())()
                    );
                }
                else if (tensorFields_.found(iter().keyword()))
                {
                    writeEntry
                    (
                        os,
                        iter().keyword(),
                        *tensorFields_.find(iter().keyword())()
                    );
                }
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
