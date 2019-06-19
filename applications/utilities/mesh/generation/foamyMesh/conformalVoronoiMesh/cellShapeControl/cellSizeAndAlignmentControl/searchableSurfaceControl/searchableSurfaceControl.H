/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::searchableSurfaceControl

Description

SourceFiles
    searchableSurfaceControl.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceControl_H
#define searchableSurfaceControl_H

#include "cellSizeFunction.H"
#include "triad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class surfaceControl Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceControl
:
    public cellSizeAndAlignmentControl
{
    // Private Data

        //- Name of the surface
        const word surfaceName_;

        //- Reference to the searchableSurface object holding the geometry data
        const searchableSurface& searchableSurface_;

        const conformationSurfaces& geometryToConformTo_;

        PtrList<cellSizeFunction> cellSizeFunctions_;

        labelList regionToCellSizeFunctions_;

        label maxPriority_;


public:

    //- Runtime type information
    TypeName("searchableSurfaceControl");


    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
        searchableSurfaceControl
        (
            const Time& runTime,
            const word& name,
            const dictionary& controlFunctionDict,
            const conformationSurfaces& geometryToConformTo,
            const scalar& defaultCellSize
        );

        //- Disallow default bitwise copy construction
        searchableSurfaceControl(const searchableSurfaceControl&) = delete;


    //- Destructor
    ~searchableSurfaceControl();


    // Member Functions

        virtual void cellSizeFunctionVertices
        (
            DynamicList<Foam::point>& pts,
            DynamicList<scalar>& sizes
        ) const;

        virtual void initialVertices
        (
            pointField& pts,
            scalarField& sizes,
            triadField& alignments
        ) const;

        const PtrList<cellSizeFunction>& sizeFunctions() const
        {
            return cellSizeFunctions_;
        }

        virtual label maxPriority() const
        {
            return maxPriority_;
        }

        bool cellSize
        (
            const Foam::point& pt,
            scalar& cellSize,
            label& priority
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaceControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
