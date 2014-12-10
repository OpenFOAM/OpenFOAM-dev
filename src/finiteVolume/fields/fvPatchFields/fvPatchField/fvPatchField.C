/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "IOobject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const word& patchType
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(patchType)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (&iF && iF.size())
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }
    this->map(ptf, mapper);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else if (!valueRequired)
    {
        fvPatchField<Type>::operator=(pTraits<Type>::zero);
    }
    else
    {
        FatalIOErrorIn
        (
            "fvPatchField<Type>::fvPatchField"
            "("
            "const fvPatch& p,"
            "const DimensionedField<Type, volMesh>& iF,"
            "const dictionary& dict,"
            "const bool valueRequired"
            ")",
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::fvPatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh();
}


template<class Type>
void Foam::fvPatchField<Type>::check(const fvPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorIn("PatchField<Type>::check(const fvPatchField<Type>&)")
            << "different patches for fvPatchField<Type>s"
            << abort(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fvPatchField<Type>::snGrad() const
{
    return patch_.deltaCoeffs()*(*this - patchInternalField());
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::fvPatchField<Type>::patchInternalField() const
{
    return patch_.patchInternalField(internalField_);
}


template<class Type>
void Foam::fvPatchField<Type>::patchInternalField(Field<Type>& pif) const
{
    patch_.patchInternalField(internalField_, pif);
}


template<class Type>
void Foam::fvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    Field<Type>& f = *this;

    if (!this->size())
    {
        f.setSize(mapper.size());
        if (f.size())
        {
            f = this->patchInternalField();
        }
    }
    else
    {
        // Map all faces provided with mapping data
        Field<Type>::autoMap(mapper);

        // For unmapped faces set to internal field value (zero-gradient)
        if
        (
            mapper.direct()
         && &mapper.directAddressing()
         && mapper.directAddressing().size()
        )
        {
            Field<Type> pif(this->patchInternalField());

            const labelList& mapAddressing = mapper.directAddressing();

            forAll(mapAddressing, i)
            {
                if (mapAddressing[i] < 0)
                {
                    f[i] = pif[i];
                }
            }
        }
        else if (!mapper.direct() && mapper.addressing().size())
        {
            Field<Type> pif(this->patchInternalField());

            const labelListList& mapAddressing = mapper.addressing();

            forAll(mapAddressing, i)
            {
                const labelList& localAddrs = mapAddressing[i];

                if (!localAddrs.size())
                {
                    f[i] = pif[i];
                }
            }
        }
    }
}


template<class Type>
void Foam::fvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvPatchField<Type>::updateCoeffs()
{
    updated_ = true;
}


template<class Type>
void Foam::fvPatchField<Type>::updateCoeffs(const scalarField& weights)
{
    if (!updated_)
    {
        updateCoeffs();

        Field<Type>& fld = *this;
        fld *= weights;

        updated_ = true;
    }
}


template<class Type>
void Foam::fvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated_)
    {
        updateCoeffs();
    }

    updated_ = false;
    manipulatedMatrix_ = false;
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix(fvMatrix<Type>& matrix)
{
    manipulatedMatrix_ = true;
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const scalarField& weights
)
{
    manipulatedMatrix_ = true;
}


template<class Type>
void Foam::fvPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;

    if (patchType_.size())
    {
        os.writeKeyword("patchType") << patchType_
            << token::END_STATEMENT << nl;
    }
}


template<class Type>
template<class EntryType>
void Foam::fvPatchField<Type>::writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const EntryType& value1,
    const EntryType& value2
) const
{
    if (value1 != value2)
    {
        os.writeKeyword(entryName) << value2 << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator*=(const fvPatchField<scalar>& ptf)"
        )   << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator/=(const fvPatchField<scalar>& ptf)"
        )   << "    incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const fvPatchField<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "fvPatchFieldNew.C"

// ************************************************************************* //
