/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "IntegrationScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::IntegrationScheme<Type>::IntegrationScheme()
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

template<class Type>
Foam::IntegrationScheme<Type>::~IntegrationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::IntegrationScheme<Type>::delta
(
    const Type& phi,
    const scalar dt,
    const Type& Alpha,
    const scalar Beta
) const
{
    return delta(phi, dt, Beta, Alpha, Beta);
}


template<class Type>
Type Foam::IntegrationScheme<Type>::delta
(
    const Type& phi,
    const scalar dt,
    const scalar Beta,
    const Type& alpha,
    const scalar beta
) const
{
    return (alpha - beta*phi)*factor(dt, Beta);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IntegrationSchemeNew.C"

// ************************************************************************* //
