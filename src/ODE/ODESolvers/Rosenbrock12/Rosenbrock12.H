/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::Rosenbrock12

Description
    L-stable embedded Rosenbrock ODE solver of order (1)2.

    References:
    \verbatim
        Verwer, J. G., Spee, E. J., Blom, J. G., & Hundsdorfer, W. (1999).
        A second-order Rosenbrock method applied to
        photochemical dispersion problems.
        SIAM Journal on Scientific Computing, 20(4), 1456-1480.
    \endverbatim

SourceFiles
    Rosenbrock12.C

\*---------------------------------------------------------------------------*/

#ifndef Rosenbrock12_H
#define Rosenbrock12_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Rosenbrock12 Declaration
\*---------------------------------------------------------------------------*/

class Rosenbrock12
:
    public ODESolver,
    public adaptiveSolver
{
    // Private Data

        mutable scalarField k1_;
        mutable scalarField k2_;
        mutable scalarField err_;
        mutable scalarField dydx_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalarSquareMatrix a_;
        mutable labelList pivotIndices_;

        static const scalar
            a21,
            c21,
            b1, b2,
            gamma,
            c2,
            e1, e2,
            d1, d2;


public:

    //- Runtime type information
    TypeName("Rosenbrock12");


    // Constructors

        //- Construct from ODESystem
        Rosenbrock12(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~Rosenbrock12()
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
            const label li,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const;

        //- Solve the ODE system and the update the state
        virtual void solve
        (
            scalar& x,
            scalarField& y,
            const label li,
            scalar& dxTry
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
