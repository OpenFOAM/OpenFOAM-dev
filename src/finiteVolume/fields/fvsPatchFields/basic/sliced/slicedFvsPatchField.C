/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "slicedFvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::slicedFvsPatchField<Type>::slicedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& completeField
)
:
    fvsPatchField<Type>(p, iF, Field<Type>())
{
    // Set the fvsPatchField to a slice of the given complete field
    UList<Type>::shallowCopy(p.patchSlice(completeField));
}


template<class Type>
Foam::slicedFvsPatchField<Type>::slicedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvsPatchField<Type>& pf
)
:
    fvsPatchField<Type>(p, iF, Field<Type>())
{
    // Set the fvsPatchField values to the given fvsPatchField
    UList<Type>::shallowCopy(pf);
}


template<class Type>
Foam::slicedFvsPatchField<Type>::slicedFvsPatchField
(
    const slicedFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf.patch(), iF, Field<Type>())
{
    // Transfer the slice from the argument
    UList<Type>::shallowCopy(ptf);
}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::slicedFvsPatchField<Type>::clone() const
{
    return tmp<fvsPatchField<Type>>
    (
        new slicedFvsPatchField<Type>(*this)
    );
}


template<class Type>
Foam::slicedFvsPatchField<Type>::slicedFvsPatchField
(
    const slicedFvsPatchField<Type>& ptf
)
:
    fvsPatchField<Type>
    (
        ptf.patch(),
        ptf.internalField(),
        Field<Type>()
    )
{
    // Transfer the slice from the argument
    UList<Type>::shallowCopy(ptf);
}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::slicedFvsPatchField<Type>::clone
(
    const DimensionedField<Type, surfaceMesh>& iF
) const
{
    return tmp<fvsPatchField<Type>>
    (
        new slicedFvsPatchField<Type>(*this, iF)
    );
}


template<class Type>
Foam::slicedFvsPatchField<Type>::~slicedFvsPatchField()
{
    // Set the fvsPatchField storage pointer to nullptr before its destruction
    // to protect the field it a slice of.
    UList<Type>::shallowCopy(UList<Type>(nullptr, 0));
}


// ************************************************************************* //
