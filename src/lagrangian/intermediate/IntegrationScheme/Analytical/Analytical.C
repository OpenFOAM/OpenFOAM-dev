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

#include "Analytical.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Analytical<Type>::Analytical
(
    const word& phiName,
    const dictionary& dict
)
:
    IntegrationScheme<Type>(phiName, dict)
{}


template<class Type>
Foam::Analytical<Type>::Analytical(const Analytical& is)
:
    IntegrationScheme<Type>(is)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Analytical<Type>::~Analytical()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::IntegrationScheme<Type>::integrationResult
Foam::Analytical<Type>::integrate
(
    const Type& phi,
    const scalar dt,
    const Type& alphaBeta,
    const scalar beta
) const
{
    typename IntegrationScheme<Type>::integrationResult retValue;

    const scalar expTerm = exp(min(50, -beta*dt));

    if (beta > ROOTVSMALL)
    {
        const Type alpha = alphaBeta/beta;
        retValue.average() = alpha + (phi - alpha)*(1 - expTerm)/(beta*dt);
        retValue.value() =  alpha + (phi - alpha)*expTerm;
    }
    else
    {
        retValue.average() = phi;
        retValue.value() =  phi;
    }


    return retValue;
}


// ************************************************************************* //
