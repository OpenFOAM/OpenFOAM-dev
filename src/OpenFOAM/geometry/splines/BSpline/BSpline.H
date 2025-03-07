/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::BSpline

Description
    An implementation of B-splines. This is a piecewise cubic spline that
    intersects its given end points, but does not intersect and of its
    intermediate control points.

    In matrix form, the \e local interpolation on the interval t=[0..1] is
    described as follows:

    \verbatim
    P(t) = 1/6 * [ t^3 t^2 t 1 ] * [ -1  3 -3  1 ] * [ P-1 ]
                                   [  3 -6  3  0 ]   [ P0 ]
                                   [ -3  0  3  0 ]   [ P1 ]
                                   [  1  4  1  0 ]   [ P2 ]
    \endverbatim

    Where P-1 and P2 represent the neighbouring points or the extrapolated
    end points. Reflection is used to automatically create the extrapolated end
    points.

See also
    CatmullRomSpline

SourceFiles
    BSpline.C

\*---------------------------------------------------------------------------*/

#ifndef BSpline_H
#define BSpline_H

#include "Spline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BSpline Declaration
\*---------------------------------------------------------------------------*/

class BSpline
:
    public Spline<BSpline>
{
public:

    // Constructors

        //- Inherit constructors
        using Spline::Spline;


    // Member Functions

        //- Return the point position corresponding to the segment parameter
        point position(const label segmenti, const scalar segmentLambda) const;

        //- Inherit the base class' position methods
        using Spline<BSpline>::position;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
