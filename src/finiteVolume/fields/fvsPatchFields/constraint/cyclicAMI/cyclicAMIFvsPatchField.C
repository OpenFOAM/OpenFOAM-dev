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

#include "cyclicAMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const cyclicAMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not cyclicAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const cyclicAMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvsPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


// ************************************************************************* //
