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
    Foam::SIBS

Description
    A semi-implicit mid-point solver for stiff systems of ordinary differential
    equations.

    Reference:
    \verbatim
        Bader, G., & Deuflhard, P. (1983).
        A semi-implicit mid-point rule for stiff systems
        of ordinary differential equations.
        Numerische Mathematik, 41(3), 373-398.
    \endverbatim

SourceFiles
    SIMPR.C
    polyExtrapolate.C

\*---------------------------------------------------------------------------*/

#ifndef SIBS_H
#define SIBS_H

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class SIBS Declaration
\*---------------------------------------------------------------------------*/

class SIBS
:
    public ODESolver
{
    // Private Data

        static const label kMaxX_ = 7, iMaxX_ = kMaxX_ + 1;
        static const label nSeq_[iMaxX_];

        static const scalar safe1, safe2, redMax, redMin, scaleMX;

        mutable scalarField a_;
        mutable scalarSquareMatrix alpha_;
        mutable scalarRectangularMatrix d_p_;
        mutable scalarField x_p_;
        mutable scalarField err_;

        mutable scalarField yTemp_;
        mutable scalarField ySeq_;
        mutable scalarField yErr_;
        mutable scalarField dydx0_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;

        mutable label first_, kMax_, kOpt_;
        mutable scalar epsOld_, xNew_;


    // Private Member Functions

        void SIMPR
        (
            const scalar xStart,
            const scalarField& y,
            const label li,
            const scalarField& dydx,
            const scalarField& dfdx,
            const scalarSquareMatrix& dfdy,
            const scalar deltaX,
            const label nSteps,
            scalarField& yEnd
        ) const;

        void polyExtrapolate
        (
            const label iest,
            const scalar xest,
            const scalarField& yest,
            scalarField& yz,
            scalarField& dy,
            scalarField& x_p,
            scalarRectangularMatrix& d_p
        ) const;


public:

    //- Runtime type information
    TypeName("SIBS");


    // Constructors

        //- Construct from ODE system
        SIBS(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~SIBS()
    {}


    // Member Functions

        //- Resize the ODE solver
        virtual bool resize();

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
