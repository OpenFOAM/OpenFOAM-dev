/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::fv::gradientLimiters::Venkatakrishnan

Description
    Venkatakrishnan gradient limiter

    to be used with the Foam::fv::cellLimitedGrad limited gradient.

    Reference:
    \verbatim
        Venkatakrishnan, V. (1993).
        On the accuracy of limiters and convergence to steady state solutions.
        In 31st Aerospace Sciences Meeting (p. 880).
    \endverbatim

    Example:
    \verbatim
    gradSchemes
    {
        default Gauss linear;
        limited cellLimited<Venkatakrishnan> Gauss linear 1;
    }
    \endverbatim

    Note: this limiter formally allows the limiter function to slightly exceed 1
    which is generally not a good idea and can cause unboundedness.  To avoid
    this problem the limiter function is clipped to 1 which is likely to be
    better behaved but is then not differentiable and so no longer conforms to
    the basic principles of this kind of limiter function.  All these issues are
    resolved in a consistent and differentiable manner by the
    Foam::fv::gradientLimiters::cubic limiter which is recommended in
    preference to the Foam::fv::gradientLimiters::Venkatakrishnan limiter.

See also
    Foam::fv::cellLimitedGrad
    Foam::fv::gradientLimiters::cubic

\*---------------------------------------------------------------------------*/

#ifndef VenkatakrishnanGradientLimiter_H
#define VenkatakrishnanGradientLimiter_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

namespace gradientLimiters
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Venkatakrishnan
{

public:

    // Constructors

        Venkatakrishnan(Istream&)
        {}


    // Member Functions

        inline scalar limiter(const scalar r) const
        {
            return (sqr(r) + 2*r)/(sqr(r) + r + 2);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace gradientLimiters

} // End namespace fv

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
