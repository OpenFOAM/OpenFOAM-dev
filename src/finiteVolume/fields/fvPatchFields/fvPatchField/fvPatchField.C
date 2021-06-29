/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    manipulatedMatrix_(false)
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
    manipulatedMatrix_(false)
{}


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
    manipulatedMatrix_(false)
{
    if (valueRequired)
    {
        if (dict.found("value"))
        {
            Field<Type>::operator=
            (
                Field<Type>("value", dict, p.size())
            );
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Essential entry 'value' missing"
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false)
{
    if (mappingRequired)
    {
        // For unmapped faces set to internal field value (zero-gradient)
        if (notNull(iF) && mapper.hasUnmapped())
        {
            fvPatchField<Type>::operator=(this->patchInternalField());
        }
        mapper(*this, ptf);
    }
}


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
    manipulatedMatrix_(false)
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
        FatalErrorInFunction
            << "different patches for fvPatchField<Type>s"
            << abort(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatchField<Type>::snGrad() const
{
    return patch_.deltaCoeffs()*(*this - patchInternalField());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
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
    mapper(*this, *this);
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
void Foam::fvPatchField<Type>::updateWeightedCoeffs(const scalarField& weights)
{
    // Default behaviour ignores the weights
    if (!updated_)
    {
        updateCoeffs();

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
    writeEntry(os, "type", type());

    if (overridesConstraint())
    {
        writeEntry(os, "patchType", patch().type());
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
        FatalErrorInFunction
            << "incompatible patches for patch fields"
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
        FatalErrorInFunction
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

#include "fvPatchFieldNew.C"

// ************************************************************************* //
