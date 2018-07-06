/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::fv::gradientLimiters::cubic

Description
    Cubic gradient limiter

    to be used with the Foam::fv::cellLimitedGrad limited gradient.  This
    limiter function is similar to Foam::fv::gradientLimiters::Venkatakrishnan
    but is a fit to obey the value and gradient constraints and avoids the
    problem of the limiter exceeding 1 present in the Venkatakrishnan function.

    The transition point at which the limiter function reaches 1 is an input
    parameter and should be set to a value between 1 and 2 although values
    larger than 2 are physical but likely to significantly reduce the accuracy
    of the scheme.

    Reference:
    \verbatim
        Michalak, K., & Ollivier-Gooch, C. (2008).
        Limiters for unstructured higher-order accurate solutions
        of the Euler equations.
        In 46th AIAA Aerospace Sciences Meeting and Exhibit (p. 776).
    \endverbatim

    Example:
    \verbatim
    gradSchemes
    {
        default Gauss linear;
        limited cellLimited<cubic> 1.5 Gauss linear 1;
    }
    \endverbatim

See also
    Foam::fv::cellLimitedGrad
    Foam::fv::gradientLimiters::Venkatakrishnan

\*---------------------------------------------------------------------------*/

#ifndef cubicGradientLimiter_H
#define cubicGradientLimiter_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

namespace gradientLimiters
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class cubic
{
    // Private Data

        //- Limiter transition point at which the limiter function -> 1
        //  Must be > 1
        const scalar rt_;

        //- Coefficient of the r^3 term (evaluated from rt)
        const scalar a_;

        // - Coefficient of the r^2 term  (evaluated from rt)
        const scalar b_;


public:

    // Constructors

        cubic(Istream& schemeData)
        :
            rt_(readScalar(schemeData)),
            a_(2.0/sqr(rt_) - 2.0/pow3(rt_)),
            b_(-(3.0/2.0)*a_*rt_)
        {
            if (rt_ < 1)
            {
                FatalIOErrorInFunction
                (
                    schemeData
                )   << "coefficient = " << rt_
                    << " should be > 1"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        inline scalar limiter(const scalar r) const
        {
            if (r < rt_)
            {
                return ((a_*r + b_)*r + 1)*r;
            }
            else
            {
                return 1;
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace gradientLimiters

} // End namespace fv

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
