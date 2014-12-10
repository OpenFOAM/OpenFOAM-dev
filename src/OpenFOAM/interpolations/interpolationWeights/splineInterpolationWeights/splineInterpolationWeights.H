/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::splineInterpolationWeights

Description
    Catmull-Rom spline interpolation.

SourceFiles
    splineInterpolationWeights.C

\*---------------------------------------------------------------------------*/

#ifndef splineInterpolationWeights_H
#define splineInterpolationWeights_H

#include "interpolationWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class splineInterpolationWeights Declaration
\*---------------------------------------------------------------------------*/

class splineInterpolationWeights
:
    public interpolationWeights
{

private:

    // Private data

        //- Cached index in samples from previous invocation
        mutable label index_;

public:

    //- Runtime type information
    TypeName("spline");

    // Constructors

        //- Construct from components. By default make sure samples are
        //  equidistant.
        splineInterpolationWeights
        (
            const scalarField& samples,
            const bool checkEqualDistance = true
        );


    //- Destructor
    virtual ~splineInterpolationWeights()
    {}


    // Member Functions

        //- Calculate weights and indices to calculate t from samples.
        //  Returns true if indices changed.
        virtual bool valueWeights
        (
            const scalar t,
            labelList& indices,
            scalarField& weights
        ) const;

        //- Calculate weights and indices to calculate integrand of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool integrationWeights
        (
            const scalar t1,
            const scalar t2,
            labelList& indices,
            scalarField& weights
        ) const
        {
            notImplemented
            (
                "splineInterpolationWeights::integrationWeights(..)"
            );
            return false;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
