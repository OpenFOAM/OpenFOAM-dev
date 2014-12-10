/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "cyclicACMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField
(
    const cyclicACMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{
    if (!isA<cyclicACMIFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField\n"
            "("
                "const cyclicACMIFvsPatchField<Type>&, "
                "const fvPatch&, "
                "const DimensionedField<Type, surfaceMesh>&, "
                "const fvPatchFieldMapper&"
            ")"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{
    if (!isA<cyclicACMIFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField"
            "("
                "const fvPatch&, "
                "const Field<Type>&, "
                "const dictionary&"
            ")",
            dict
        )   << "patch " << this->patch().index() << " not cyclicACMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField
(
    const cyclicACMIFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_)
{}


template<class Type>
Foam::cyclicACMIFvsPatchField<Type>::cyclicACMIFvsPatchField
(
    const cyclicACMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicACMIFvsPatchField<Type>::coupled() const
{
    if
    (
        Pstream::parRun()
     || (
            this->cyclicACMIPatch_.size()
         && this->cyclicACMIPatch_.cyclicACMIPatch().neighbPatch().size()
        )
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
