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
    Foam::Spline

Description
    Base class for splines

SourceFiles
    Spline.C
    SplineI.H

\*---------------------------------------------------------------------------*/

#ifndef Spline_H
#define Spline_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Spline Declaration
\*---------------------------------------------------------------------------*/

template<class Derived>
class Spline
{
    // Private Data

        //- The end points and control points
        const pointField knots_;

        //- Tolerance to which to improve the spacing
        const scalar tol_;

        //- Number of iterations over which to improve the spacing
        const label nIter_;


protected:

    // Protected Member Functions

        //- Access the end points and control points
        inline const pointField& knots() const;


public:

    // Constructors

        //- Construct from components
        Spline(const pointField& knots, const scalar tol, const label nIter);

        //- Construct from knots with default iteration settings
        Spline(const pointField& knots);

        //- Disallow default bitwise copy construction
        Spline(const Spline&) = delete;


    // Member Functions

        //- Return the number of segments
        inline label nSegments() const;

        //- Return the start point
        inline const point& start() const;

        //- Return the end point
        inline const point& end() const;

        //- Return the point position corresponding to the curve parameter
        point position(const scalar lambda) const;

        //- Return the point positions corresponding to the curve parameters
        tmp<pointField> position(const scalarList& lambdas) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Spline&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SplineI.H"

#ifdef NoRepository
    #include "Spline.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
