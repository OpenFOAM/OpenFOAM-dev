/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "internalFvFieldSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::internalFvFieldSource<Type>::internalFvFieldSource
(
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvFieldSource<Type>(iF, dict)
{}


template<class Type>
Foam::internalFvFieldSource<Type>::internalFvFieldSource
(
    const internalFvFieldSource<Type>& field,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvFieldSource<Type>(field, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::internalFvFieldSource<Type>::~internalFvFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalFvFieldSource<Type>::sourceValue(const fvSource& source) const
{
    // This value doesn't matter in principle, as this condition takes 100% of
    // its value from the internal field. However, this value does, at least,
    // need to be physical to prevent things like thermo evaluations from
    // failing. So, just take the internal values.
    return
        tmp<Field<Type>>
        (
            new Field<Type>(this->internalField(), source.cells())
        );
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::internalFvFieldSource<Type>::internalCoeff(const fvSource& source) const
{
    return tmp<scalarField>(new scalarField(source.nCells(), scalar(1)));
}


// ************************************************************************* //
