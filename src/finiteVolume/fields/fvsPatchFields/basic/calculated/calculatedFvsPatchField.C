/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "calculatedFvsPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
const word& fvsPatchField<Type>::calculatedType()
{
    return calculatedFvsPatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(p, iF)
{}


template<class Type>
calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf
)
:
    fvsPatchField<Type>(ptf)
{}


template<class Type>
calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf, iF)
{}


template<class Type>
tmp<fvsPatchField<Type> > fvsPatchField<Type>::NewCalculatedType
(
    const fvPatch& p
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != patchConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()
        (
            p,
            DimensionedField<Type, surfaceMesh>::null()
        );
    }
    else
    {
        return tmp<fvsPatchField<Type> >
        (
            new calculatedFvsPatchField<Type>
            (
                p,
                DimensionedField<Type, surfaceMesh>::null()
            )
        );
    }
}


template<class Type>
template<class Type2>
tmp<fvsPatchField<Type> > fvsPatchField<Type>::NewCalculatedType
(
    const fvsPatchField<Type2>& pf
)
{
    return NewCalculatedType(pf.patch());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
