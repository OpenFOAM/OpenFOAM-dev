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

#include "calculatedPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
const word& pointPatchField<Type>::calculatedType()
{
    return calculatedPointPatchField<Type>::typeName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(p, iF)
{}


template<class Type>
calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchField<Type>(p, iF, dict)
{}


template<class Type>
calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const calculatedPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    pointPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const calculatedPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
autoPtr<pointPatchField<Type> >
pointPatchField<Type>::NewCalculatedType
(
    const pointPatchField<Type2>& pf
)
{
    typename pointPatchConstructorTable::iterator patchTypeCstrIter =
        pointPatchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != pointPatchConstructorTablePtr_->end())
    {
        return autoPtr<pointPatchField<Type> >
        (
            patchTypeCstrIter()
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
    else
    {
        return autoPtr<pointPatchField<Type> >
        (
            new calculatedPointPatchField<Type>
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
