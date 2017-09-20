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

Class
    Foam::Analytical

Description
    Analytical integration scheme

    \f[
        \Delta \phi = (A - B \phi^n) \frac{1}{B} (1 - e^{- B \Delta t})
    \f]

\*---------------------------------------------------------------------------*/

#ifndef Analytical_H
#define Analytical_H

#include "IntegrationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace integrationSchemes
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

        //- Construct
        Analytical();

        //- Construct and return clone
        virtual autoPtr<IntegrationScheme<Type>> clone() const
        {
            return autoPtr<IntegrationScheme<Type>>
            (
                new Analytical<Type>(*this)
            );
        }


    //- Destructor
    virtual ~Analytical();


    // Member Functions

        //- Return the integration factor
        virtual scalar factor(const scalar dt, const scalar Beta) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace integrationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Analytical.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
