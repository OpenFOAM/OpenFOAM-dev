/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::linearDistance

Description

SourceFiles
    linearDistance.C

\*---------------------------------------------------------------------------*/

#ifndef linearDistance_H
#define linearDistance_H

#include "cellSizeFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearDistance Declaration
\*---------------------------------------------------------------------------*/

class linearDistance
:
    public cellSizeFunction
{

private:

    // Private data

        //- cell size at distance_ from the surface
        scalar distanceCellSize_;

        //- distance from the surface to control over
        scalar distance_;

        //- distance squared
        scalar distanceSqr_;


    // Private Member Functions

        //- Calculate the cell size as a function of the given distance
        scalar sizeFunction(const point& pt, scalar d, label index) const;


public:

    //- Runtime type information
    TypeName("linearDistance");

    // Constructors

        //- Construct from components
        linearDistance
        (
            const dictionary& initialPointsDict,
            const searchableSurface& surface,
            const scalar& defaultCellSize,
            const labelList regionIndices
        );


    //- Destructor
    virtual ~linearDistance()
    {}


    // Member Functions


        virtual bool sizeLocations
        (
            const pointIndexHit& hitPt,
            const vector& n,
            pointField& shapePts,
            scalarField& shapeSizes
        ) const;

        //- Modify scalar argument to the cell size specified by function.
        //  Return a boolean specifying if the function was used, i.e. false if
        //  the point was not in range of the surface for a spatially varying
        //  size.
        virtual bool cellSize
        (
            const point& pt,
            scalar& size
        ) const;

        //- Adapt local cell size. Return true if anything changed.
        virtual bool setCellSize
        (
            const pointField& pts
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
