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
    Foam::cellShapeControl

Description

SourceFiles
    cellShapeControlI.H
    cellShapeControl.C

\*---------------------------------------------------------------------------*/

#ifndef cellShapeControl_H
#define cellShapeControl_H

#include "dictionary.H"
#include "autoPtr.H"
#include "tensor.H"
#include "point.H"
#include "primitiveFieldsFwd.H"
#include "pointFieldFwd.H"
#include "Time.H"
#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "cellAspectRatioControl.H"
#include "cellSizeAndAlignmentControls.H"
#include "cellShapeControlMesh.H"
#include "backgroundMeshDecomposition.H"
#include "cvControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class cellShapeControl Declaration
\*---------------------------------------------------------------------------*/

class cellShapeControl
:
    public dictionary
{
    // Private data

        const Time& runTime_;

        const searchableSurfaces& allGeometry_;

        const conformationSurfaces& geometryToConformTo_;

        const scalar defaultCellSize_;

        const scalar minimumCellSize_;

        cellShapeControlMesh shapeControlMesh_;

        cellAspectRatioControl aspectRatio_;

        cellSizeAndAlignmentControls sizeAndAlignment_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellShapeControl(const cellShapeControl&);

        //- Disallow default bitwise assignment
        void operator=(const cellShapeControl&);


public:

    //- Runtime type information
    ClassName("cellShapeControl");


    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
        cellShapeControl
        (
            const Time& runTime,
            const cvControls& foamyHexMeshControls,
            const searchableSurfaces& allGeometry,
            const conformationSurfaces& geometryToConformTo
        );


    //- Destructor
    ~cellShapeControl();


    // Member Functions

        // Access

            inline const scalar& defaultCellSize() const;

            inline cellShapeControlMesh& shapeControlMesh();

            inline const cellShapeControlMesh& shapeControlMesh() const;

            inline const cellAspectRatioControl& aspectRatio() const;

            inline const cellSizeAndAlignmentControls& sizeAndAlignment() const;

            inline const scalar& minimumCellSize() const;


        // Query

            //- Return the cell size at the given location
            scalar cellSize(const point& pt) const;

            scalarField cellSize(const pointField& pts) const;

            //- Return the cell alignment at the given location
            tensor cellAlignment(const point& pt) const;

            void cellSizeAndAlignment
            (
                const point& pt,
                scalar& size,
                tensor& alignment
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellShapeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
