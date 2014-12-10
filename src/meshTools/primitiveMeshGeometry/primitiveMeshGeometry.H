/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::primitiveMeshGeometry

Description
    Updateable mesh geometry + checking routines.

SourceFiles
    primitiveMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef primitiveMeshGeometry_H
#define primitiveMeshGeometry_H

#include "pointFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class primitiveMeshGeometry Declaration
\*---------------------------------------------------------------------------*/

class primitiveMeshGeometry
{
        //- Reference to primitiveMesh.
        const primitiveMesh& mesh_;

        //- Uptodate copy of face areas
        vectorField faceAreas_;

        //- Uptodate copy of face centres
        vectorField faceCentres_;

        //- Uptodate copy of cell centres
        vectorField cellCentres_;

        //- Uptodate copy of cell volumes
        scalarField cellVolumes_;


    // Private Member Functions

        //- Update face areas and centres on selected faces.
        void updateFaceCentresAndAreas
        (
            const pointField& p,
            const labelList& changedFaces
        );

        //- Update cell volumes and centres on selected cells. Requires
        //  cells and faces to be consistent set.
        void updateCellCentresAndVols
        (
            const labelList& changedCells,
            const labelList& changedFaces
        );

public:

    ClassName("primitiveMeshGeometry");

    // Constructors

        //- Construct from mesh
        primitiveMeshGeometry(const primitiveMesh&);


    // Member Functions

        // Access

            const primitiveMesh& mesh() const
            {
                return mesh_;
            }

            const vectorField& faceAreas() const
            {
                return faceAreas_;
            }
            const vectorField& faceCentres() const
            {
                return faceCentres_;
            }
            const vectorField& cellCentres() const
            {
                return cellCentres_;
            }
            const scalarField& cellVolumes() const
            {
                return cellVolumes_;
            }

        // Edit

            //- Take over properties from mesh
            void correct();

            //- Recalculate on selected faces. Recalculates cell properties
            //  on owner and neighbour of these cells.
            void correct
            (
                const pointField& p,
                const labelList& changedFaces
            );

            //- Helper function: get affected cells from faces
            labelList affectedCells(const labelList& changedFaces) const;


        // Checking of selected faces with supplied geometry (mesh only used for
        // topology). Parallel aware.

            static bool checkFaceDotProduct
            (
                const bool report,
                const scalar orthWarn,
                const primitiveMesh&,
                const vectorField& cellCentres,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            static bool checkFacePyramids
            (
                const bool report,
                const scalar minPyrVol,
                const primitiveMesh&,
                const vectorField& cellCentres,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet*
            );

            static bool checkFaceSkewness
            (
                const bool report,
                const scalar internalSkew,
                const scalar boundarySkew,
                const primitiveMesh& mesh,
                const vectorField& cellCentres,
                const vectorField& faceCentres,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            static bool checkFaceWeights
            (
                const bool report,
                const scalar warnWeight,
                const primitiveMesh& mesh,
                const vectorField& cellCentres,
                const vectorField& faceCentres,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            static bool checkFaceAngles
            (
                const bool report,
                const scalar maxDeg,
                const primitiveMesh& mesh,
                const vectorField& faceAreas,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //static bool checkFaceFlatness
            //(
            //    const bool report,
            //    const scalar warnFlatness,
            //    const primitiveMesh&,
            //    const vectorField& faceAreas,
            //    const vectorField& faceCentres,
            //    const pointField& p,
            //    const labelList& checkFaces,
            //    labelHashSet* setPtr
            //);

            static bool checkFaceTwist
            (
                const bool report,
                const scalar minTwist,
                const primitiveMesh&,
                const vectorField& faceAreas,
                const vectorField& faceCentres,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            static bool checkFaceArea
            (
                const bool report,
                const scalar minArea,
                const primitiveMesh&,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            static bool checkCellDeterminant
            (
                const bool report,
                const scalar minDet,
                const primitiveMesh&,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr
            );


        // Checking of selected faces with local geometry. Uses above static
        // functions. Parallel aware.

            bool checkFaceDotProduct
            (
                const bool report,
                const scalar orthWarn,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFacePyramids
            (
                const bool report,
                const scalar minPyrVol,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceSkewness
            (
                const bool report,
                const scalar internalSkew,
                const scalar boundarySkew,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceWeights
            (
                const bool report,
                const scalar warnWeight,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceAngles
            (
                const bool report,
                const scalar maxDeg,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            //bool checkFaceFlatness
            //(
            //    const bool report,
            //    const scalar warnFlatness,
            //    const pointField& p,
            //    const labelList& checkFaces,
            //    labelHashSet* setPtr
            //) const;

            bool checkFaceTwist
            (
                const bool report,
                const scalar minTwist,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceArea
            (
                const bool report,
                const scalar minArea,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkCellDeterminant
            (
                const bool report,
                const scalar warnDet,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
