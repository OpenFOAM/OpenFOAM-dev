/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::adaptiveSolver

Description

SourceFiles
    adaptiveSolver.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveSolver_H
#define adaptiveSolver_H

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adaptiveSolver Declaration
\*---------------------------------------------------------------------------*/

class adaptiveSolver
{
    // Private Data

        //- Step-size adjustment controls
        scalar safeScale_, alphaInc_, alphaDec_, minScale_, maxScale_;

        //- Cache for dydx at the initial time
        mutable scalarField dydx0_;

        //- Temprorary for the test-step solution
        mutable scalarField yTemp_;


public:

    // Constructors

        //- Construct from ODESystem
        adaptiveSolver(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~adaptiveSolver()
    {}


    // Member Functions

        //- Resize the ODE solver
        bool resize(const label n);

        //- Solve a single step dx and return the error
        virtual scalar solve
        (
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const = 0;

        //- Solve the ODE system and the update the state
        void solve
        (
            const ODESystem& ode,
            scalar& x,
            scalarField& y,
            scalar& dxTry
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
