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
    Foam::rodas34

Description
    L-stable, stiffly-accurate embedded Rosenbrock ODE solver of order (3)4.

    Reference:
    \verbatim
        Hairer, E., Nørsett, S. P., & Wanner, G. (1996).
        Solving Ordinary Differential Equations II:
        Stiff and Differential-Algebraic Problems, second edition",
        Springer-Verlag, Berlin.
    \endverbatim

SourceFiles
    rodas34.C

\*---------------------------------------------------------------------------*/

#ifndef rodas34_H
#define rodas34_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rodas34 Declaration
\*---------------------------------------------------------------------------*/

class rodas34
:
    public ODESolver,
    public adaptiveSolver
{
    // Private Data

        mutable scalarField k1_;
        mutable scalarField k2_;
        mutable scalarField k3_;
        mutable scalarField k4_;
        mutable scalarField k5_;
        mutable scalarField dy_;
        mutable scalarField err_;
        mutable scalarField dydx_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalarSquareMatrix a_;
        mutable labelList pivotIndices_;

        static const scalar
            c2, c3, c4,
            d1, d2, d3, d4,
            a21, a31, a32,
            a41, a42, a43,
            a51, a52, a53, a54,
            c21, c31, c32,
            c41, c42, c43,
            c51, c52, c53, c54,
            c61, c62, c63, c64, c65,
            gamma;

public:

    //- Runtime type information
    TypeName("rodas34");


    // Constructors

        //- Construct from ODESystem
        rodas34(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~rodas34()
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
