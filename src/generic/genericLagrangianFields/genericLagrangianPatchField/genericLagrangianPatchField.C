/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "genericLagrangianPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericLagrangianPatchField<Type>::genericLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    genericFieldBase(dict.lookup("type")),
    calculatedLagrangianPatchField<Type>(p, iIo, dict),
    dict_(dict)
{}


template<class Type>
Foam::genericLagrangianPatchField<Type>::genericLagrangianPatchField
(
    const genericLagrangianPatchField<Type>& ptf
)
:
    genericFieldBase(ptf),
    calculatedLagrangianPatchField<Type>(ptf),
    dict_(ptf.dict_)
{}


template<class Type>
Foam::genericLagrangianPatchField<Type>::genericLagrangianPatchField
(
    const genericLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    genericFieldBase(ptf),
    calculatedLagrangianPatchField<Type>(ptf, iIo),
    dict_(ptf.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericLagrangianPatchField<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", actualTypeName());

    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type")
        {
            iter().write(os);
        }
    }
}


// ************************************************************************* //
