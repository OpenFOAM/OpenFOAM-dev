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

#include "genericLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericLagrangianFieldSource<Type>::genericLagrangianFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    genericFieldBase(dict.lookup("type")),
    LagrangianFieldSource<Type>(iIo, dict),
    dict_(dict)
{}


template<class Type>
Foam::genericLagrangianFieldSource<Type>::genericLagrangianFieldSource
(
    const genericLagrangianFieldSource<Type>& stf,
    const regIOobject& iIo
)
:
    genericFieldBase(stf),
    LagrangianFieldSource<Type>(stf, iIo),
    dict_(stf.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::genericLagrangianFieldSource<Type>::~genericLagrangianFieldSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericLagrangianFieldSource<Type>::write(Ostream& os) const
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
