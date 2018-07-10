/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    // Private data

        //- Name of the surface
        const word surfaceName_;

        //- Reference to the searchableSurface object holding the geometry data
        const searchableSurface& searchableSurface_;

        const conformationSurfaces& geometryToConformTo_;

        PtrList<cellSizeFunction> cellSizeFunctions_;

        labelList regionToCellSizeFunctions_;

        label maxPriority_;


//        const conformationSurfaces& geometryToConformTo_;
//
//        //- Indices of surfaces in allGeometry that are to be conformed to
//        labelList surfaces_;
//
//        //- A list of all of the cellSizeFunction objects
//        PtrList<cellSizeFunction> cellSizeFunctions_;
//
//        autoPtr<triangulatedMesh> triangulatedMesh_;
//
//
//    // Private Member Functions
//
//        //-
//        tensor requiredAlignment
//        (
//            const point& pt,
//            const vectorField& ptNormals
//        ) const;

        //- Disallow default bitwise copy construct
        searchableSurfaceControl(const searchableSurfaceControl&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaceControl&);


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

    //- Destructor
    ~searchableSurfaceControl();


    // Member Functions

        // Access

//            //- Return reference to the searchableSurfaces object containing
//            //  all of the geometry
//            inline const searchableSurfaces& geometry() const;
//
//            //- Return the surface indices
//            inline const labelList& surfaces() const;
//
//
//        // Query
//
//            //- Return the cell size at the given location
//            virtual scalar cellSize(const point& pt) const;
//
//            //- Return the cell alignment at the given location
//            virtual tensor cellAlignment(const point& pt) const;
//
//            virtual void cellSizeAndAlignment
//            (
//                const point& pt,
//                scalar& size,
//                tensor& alignment
//            ) const;

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

        // Edit

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
