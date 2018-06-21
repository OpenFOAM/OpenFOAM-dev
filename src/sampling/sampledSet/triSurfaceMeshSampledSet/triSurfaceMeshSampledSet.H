/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::triSurfaceMeshSampledSet

Description
    Samples from all the points of a triSurfaceMesh. Surface files are read
    from constant/triSurface.

Usage
    \table
        Property    | Description                            | Req'd? | Default
        surface     | The surface file name                  | yes    |
        axis        | The coordinate axis that is written    | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        triSurfaceMesh;
        surface     "surface.stl";
        axis        x;
    }
    \endverbatim

SourceFiles
    triSurfaceMeshSampledSet.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMeshSampledSet_H
#define triSurfaceMeshSampledSet_H

#include "sampledSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                   Class triSurfaceMeshSampledSet Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMeshSampledSet
:
    public sampledSet
{
    // Private data

        //- Name of triSurfaceMesh
        const word surface_;

        //- Sampling points
        const List<point> points_;


    // Private Member Functions

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("triSurfaceMesh");


    // Constructors

        //- Construct from dictionary
        triSurfaceMeshSampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~triSurfaceMeshSampledSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
