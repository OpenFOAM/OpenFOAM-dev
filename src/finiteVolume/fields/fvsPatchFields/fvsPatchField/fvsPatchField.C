/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(ptf, mapper),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF)
{
    if (dict.found("value"))
    {
        fvsPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "fvsPatchField<Type>::fvsPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, surfaceMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "essential value entry not provided"
            << exit(FatalIOError);
    }
}


template<class Type>
fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_)
{}


template<class Type>
fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const objectRegistry& fvsPatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh();
}


template<class Type>
void fvsPatchField<Type>::check(const fvsPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorIn("PatchField<Type>::check(const fvsPatchField<Type>&)")
            << "different patches for fvsPatchField<Type>s"
            << abort(FatalError);
    }
}


// Map from self
template<class Type>
void fvsPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


// Reverse-map the given fvsPatchField onto this fvsPatchField
template<class Type>
void fvsPatchField<Type>::rmap
(
    const fvsPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


// Write
template<class Type>
void fvsPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void fvsPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void fvsPatchField<Type>::operator=
(
    const fvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void fvsPatchField<Type>::operator+=
(
    const fvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void fvsPatchField<Type>::operator-=
(
    const fvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void fvsPatchField<Type>::operator*=
(
    const fvsPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator*=(const fvsPatchField<scalar>& ptf)"
        )   << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void fvsPatchField<Type>::operator/=
(
    const fvsPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator/=(const fvsPatchField<scalar>& ptf)"
        )   << "    incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void fvsPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void fvsPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void fvsPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void fvsPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void fvsPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void fvsPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void fvsPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void fvsPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void fvsPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void fvsPatchField<Type>::operator==
(
    const fvsPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void fvsPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void fvsPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const fvsPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const fvsPatchField<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "fvsPatchFieldNew.C"

// ************************************************************************* //
