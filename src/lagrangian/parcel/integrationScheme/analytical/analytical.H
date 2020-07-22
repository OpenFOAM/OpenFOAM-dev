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

Class
    Foam::integrationSchemes::analytical

Description
    Analytical integration scheme

    \f[
        \Delta \phi = (A - B \phi^n) \frac{1}{B} (1 - e^{- B \Delta t})
    \f]

\*---------------------------------------------------------------------------*/

#ifndef analytical_H
#define analytical_H

#include "integrationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace integrationSchemes
{

/*---------------------------------------------------------------------------*\
                         Class analytical Declaration
\*---------------------------------------------------------------------------*/

class analytical
:
    public integrationScheme
{
public:

    //- Runtime type information
    TypeName("analytical");


    // Constructors

        //- Construct
        analytical();

        //- Construct and return clone
        virtual autoPtr<integrationScheme> clone() const
        {
            return autoPtr<integrationScheme>(new analytical(*this));
        }


    //- Destructor
    virtual ~analytical();


    // Member Functions

        //- Return the integration effective time step
        virtual scalar dtEff(const scalar dt, const scalar Beta) const;

        //- Return the integral of the effective time step
        virtual scalar sumDtEff(const scalar dt, const scalar Beta) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace integrationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
