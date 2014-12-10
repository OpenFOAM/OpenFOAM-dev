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
    Foam::autoDensity

Description
    Choose random points inside the domain and place them with a probability
    proportional to the target density of points.

SourceFiles
    autoDensity.C

\*---------------------------------------------------------------------------*/

#ifndef autoDensity_H
#define autoDensity_H

#include "initialPointsMethod.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class autoDensity Declaration
\*---------------------------------------------------------------------------*/

class autoDensity
:
    public initialPointsMethod
{

private:

    // Private data

        //- Trial points attempted to be placed in all boxes
        mutable label globalTrialPoints_;

        //- Smallest minimum cell size allowed, i.e. to avoid high initial
        //  population of areas of small size
        scalar minCellSizeLimit_;

        //- Minimum normal level of recursion, can be more if a high density
        //  ratio is detected
        label minLevels_;

        //- Maximum allowed ratio of cell size in a box
        scalar maxSizeRatio_;

        //- How fine should the initial sample of the volume a box be to
        //  investigate its cell sizes and volume fraction
        label volRes_;

        //- How fine should the initial sample of the surface of a box be to
        //  investigate if it is near to a the geometry.
        label surfRes_;


    // Private Member Functions

        //- Check if the given box overlaps the geometry or, in parallel, the
        //  backgroundMeshDecomposition
        bool combinedOverlaps(const treeBoundBox& box) const;

        //- Check if the given point is inside the geometry and, in parallel,
        //  the backgroundMeshDecomposition
        bool combinedInside(const point& p) const;

        //- Check if the given points are wellInside the geometry and, in
        //  parallel, inside the backgroundMeshDecomposition
        Field<bool> combinedWellInside
        (
            const pointField& pts,
            const scalarField& sizes
        ) const;

        //- Check if the given points are wellInside the geometry and, in
        //  parallel, inside the backgroundMeshDecomposition
        bool combinedWellInside
        (
            const point& p,
            scalar size
        ) const;

        //- Write boundBox as obj
        void writeOBJ
        (
            const treeBoundBox& bb,
            fileName name
        ) const;

        //- Descend into octants of the supplied bound box
        label recurseAndFill
        (
            DynamicList<Vb::Point>& initialPoints,
            const treeBoundBox& bb,
            label levelLimit,
            word recursionName
        ) const;

        //- Fill the given box, optionally filling surface overlapping boxes.
        //  Returns true if the fill is successful, false if it is to be aborted
        //  in favour of further recursion.
        bool fillBox
        (
            DynamicList<Vb::Point>& initialPoints,
            const treeBoundBox& bb,
            bool overlapping
        ) const;


public:

    //- Runtime type information
    TypeName("autoDensity");

    // Constructors

        //- Construct from components
        autoDensity
        (
            const dictionary& initialPointsDict,
            const Time& runTime,
            Random& rndGen,
            const conformationSurfaces& geometryToConformTo,
            const cellShapeControl& cellShapeControls,
            const autoPtr<backgroundMeshDecomposition>& decomposition
        );


    //- Destructor
    virtual ~autoDensity()
    {}


    // Member Functions

        //- Return the initial points for the conformalVoronoiMesh
        virtual List<Vb::Point> initialPoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
