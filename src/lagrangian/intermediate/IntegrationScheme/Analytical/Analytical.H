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

Class
    Foam::Analytical

Description
    Analytical integration

\*---------------------------------------------------------------------------*/

#ifndef Analytical_H
#define Analytical_H

#include "IntegrationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Analytical Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Analytical
:
    public IntegrationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("analytical");


    // Constructors

        //- Construct from components
        Analytical(const word& phiName, const dictionary& dict);

        //- Copy constructor
        Analytical(const Analytical& is);

        //- Construct and return clone
        virtual autoPtr<IntegrationScheme<Type> > clone() const
        {
            return autoPtr<IntegrationScheme<Type> >
            (
                new Analytical<Type>(*this)
            );
        }


    //- Destructor
    virtual ~Analytical();


    // Member Functions

        //- Perform the integration
        virtual typename IntegrationScheme<Type>::integrationResult integrate
        (
            const Type& phi,
            const scalar dt,
            const Type& alphaBeta,
            const scalar beta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Analytical.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
