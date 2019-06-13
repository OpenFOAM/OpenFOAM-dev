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
    Foam::Euler

Description
    Euler ODE solver of order (0)1.

    The method calculates the new state from:
    \f[
        y_{n+1} = y_n + \delta_x f
    \f]
    The error is estimated directly from the change in the solution,
    i.e. the difference between the 0th and 1st order solutions:
    \f[
        err_{n+1} = y_{n+1} - y_n
    \f]

SourceFiles
    Euler.C

\*---------------------------------------------------------------------------*/

#ifndef Euler_H
#define Euler_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler
:
    public ODESolver,
    public adaptiveSolver
{
    // Private Data

        mutable scalarField err_;


public:

    //- Runtime type information
    TypeName("Euler");


    // Constructors

        //- Construct from ODESystem
        Euler(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~Euler()
    {}


    // Member Functions

        //- Inherit solve from ODESolver
        using ODESolver::solve;

        //- Resize the ODE solver
        virtual bool resize();

        //- Solve a single step dx and return the error
        virtual scalar solve
        (
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const;

        //- Solve the ODE system and the update the state
        virtual void solve
        (
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
