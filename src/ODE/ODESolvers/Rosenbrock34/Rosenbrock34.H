/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::Rosenbrock34

Description
    L-stable embedded Rosenbrock ODE solver of order (3)4.

    \verbatim
        "Solving Ordinary Differential Equations II: Stiff
         and Differential-Algebraic Problems, second edition",
        Hairer, E.,
        Nørsett, S.,
        Wanner, G.,
        Springer-Verlag, Berlin. 1996.
    \endverbatim

    The default constants are from:
    \verbatim
        "Implementation of Rosenbrock Methods"
        Shampine, L. F.,
        ACM Transactions on Mathematical Software, vol. 8, 1982, pp. 93–113.
    \endverbatim
    with which the scheme is more accurate than with the L-Stable coefficients
    for small step-size but less stable for large step-size.

    The L-Stable scheme constants are provided commented-out in Rosenbrock34.C

SourceFiles
    Rosenbrock34.C

\*---------------------------------------------------------------------------*/

#ifndef Rosenbrock34_H
#define Rosenbrock34_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Rosenbrock34 Declaration
\*---------------------------------------------------------------------------*/

class Rosenbrock34
:
    public ODESolver,
    public adaptiveSolver
{
    // Private data

        mutable scalarField k1_;
        mutable scalarField k2_;
        mutable scalarField k3_;
        mutable scalarField k4_;
        mutable scalarField err_;
        mutable scalarField dydx_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalarSquareMatrix a_;
        mutable labelList pivotIndices_;

        static const scalar
            a21, a31, a32,
            c21, c31, c32,
            c41, c42, c43,
            b1, b2, b3, b4,
            e1, e2, e3, e4,
            gamma,
            c2, c3,
            d1, d2, d3, d4;


public:

    //- Runtime type information
    TypeName("Rosenbrock34");


    // Constructors

        //- Construct from ODE
        Rosenbrock34(const ODESystem& ode, const dictionary& dict);


    // Member Functions

        //- Solve a single step dx and return the error
        scalar solve
        (
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const;

        //- Solve the ODE system and the update the state
        void solve
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
