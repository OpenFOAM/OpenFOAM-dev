/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::MPLICcell

Description
    Class performs geometric matching of volume fraction and calculates surface
    interpolation of volume fraction field.

    Cut algorithms:
    - Single cell cut
    - Face-edge walk multiple cell cuts
    - Tet decomposition cell cuts

SourceFiles
    MPLICcell.C

\*---------------------------------------------------------------------------*/

#ifndef MPLICcell_H
#define MPLICcell_H

#include "MPLICface.H"
#include "MPLICcellStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MPLICcell Declaration
\*---------------------------------------------------------------------------*/

class MPLICcell
{
    // Private member data

        // Face related objects

            //- Select unweighted interpolation if true
            // velocity flux corrected if false
            const bool unweighted_;

            //- Select multi-cut if true or single-cut if false
            const bool multiCut_;

            //- Surface interpolated alpha
            DynamicList<scalar> alphaf_;

            //- Flux of alpha from point interpolated U
            DynamicList<scalar> alphaPhiU_;

            //- Face flux from point interpolated U
            DynamicList<scalar> phiU_;


        // Geometry related objects

            //- Calculated volume fraction corresponding to the cut
            scalar cutAlpha_;

            //- Cut surface area vector
            vector cutSf_;

            //- Cut normal
            vector cutNormal_;

            //- Face cutter
            MPLICface faceCutter_;

            //- Sub cell volume
            scalar subCellVolume_;

            //- Cut points for single cut
            DynamicList<point> cutPoints_;

            //- Cut edge labels for single cut
            DynamicList<label> cutEdges_;

            //- Submerged face areas
            DynamicList<vector> subFaceAreas_;

            //- Submerged face areas
            DynamicList<scalar> subFaceMagSf_;

            //- Submerged face centers
            DynamicList<vector> subFaceCentres_;

            template<class Type>
            class Vector4
            :
                public VectorSpace<Vector4<Type>, Type, 4>
            {
            public:

                // Constructors

                //- Construct null
                inline Vector4()
                {}

                inline Type a() const
                {
                    return this->v_[0];
                }

                inline Type b() const
                {
                    return this->v_[1];
                }

                inline Type c() const
                {
                    return this->v_[2];
                }

                inline Type d() const
                {
                    return this->v_[3];
                }


                inline Type& a()
                {
                    return this->v_[0];
                }

                inline Type& b()
                {
                    return this->v_[1];
                }

                inline Type& c()
                {
                    return this->v_[2];
                }

                inline Type& d()
                {
                    return this->v_[3];
                }
            };

            typedef Vector4<scalar> vector4;

            //- Four point alpha values for cubic polynomial fit
            vector4 pCubicAlphas_;

            //- Four cell alpha values for cubic polynomial fit
            vector4 cCubicAlphas_;


        // Tetrahedron storage objects

            //- Tetrahedron alpha point values
            scalarField tetPointsAlpha_;

            //- Tetrahedron velocity point values
            vectorField tetPointsU_;

            //- Tetrahedron points
            pointField tetPoints_;

            //- Tetrahedron surface area vectors
            Vector4<vector> tetSf_;

            //- Tetrahedron face centres
            Vector4<vector> tetCf_;

            //- Tetrahedron triangular faces addressing
            const FixedList<face, 4> tetFaces_;


        // Multicut addressing

            //- Is addressnig computed?
            bool addressingCalculated_;

            //- Local edge faces
            DynamicList<DynamicList<label>> localEdgeFaces_;

            //- Local face edges
            DynamicList<DynamicList<label>> localFaceEdges_;


        // Cell-point work arrays

            DynamicList<scalar> cellPointsAlpha_;

            DynamicList<scalar> pointsAlpha_;


    // Private Member Functions

        //- Match cell volume ratio to phase fraction
        //  Returns:
        //  - +1: success
        //  -  0: fail
        //  - -1: base scheme
        label calcMatchAlphaCutCell
        (
            const MPLICcellStorage& cellInfo,
            const bool tetDecom = false
        );

        //- Find in between which two point alpha values
        //  the target volume fraction lies
        void findPointAlphaBounds
        (
            const MPLICcellStorage& cellInfo,
            const bool tetDecom
        );

        //- Calculate the two interior point alpha values for the cubic fit
        void calcPointAlphaInterior
        (
            const MPLICcellStorage& cellInfo,
            const bool tetDecom
        );

        //- Solve the cubic fit
        FixedList<scalar, 4> solveVanderMatrix() const;

        //- Identifying roots of cubic equation matching target volume fraction
        void findRoots
        (
            const MPLICcellStorage& cellInfo,
            const FixedList<scalar, 4>& coeffs,
            const bool tetDecom
        );

        //- Select simple cut or tet decomposition
        scalar calcAlpha
        (
            const MPLICcellStorage& cellInfo,
            const scalar target,
            const bool tetDecom
        );

        //- Calculate current sub-cell volume
        void calcSubCellVolume();

        //- Calculate cell cut, volume and alpha
        scalar calcCutCellVolumeAlpha
        (
            const MPLICcellStorage& cellInfo,
            const scalar target
        );

        //- Calculate cell cut, volume and alpha by tet decomposition
        scalar calcTetCutCellVolumeAlpha
        (
            const MPLICcellStorage& cellInfo,
            const scalar target
        );

        //- Attempt single cut through the cell
        //  Returns:
        //  - 1: success
        //  - 0: fail
        bool singleCutCell
        (
            const MPLICcellStorage& cellInfo,
            const scalar target
        );

        //- Attempt multiple cuts through the cell
        bool multiCutCell
        (
            const MPLICcellStorage& cellInfo,
            const scalar target
        );

        //- Single cut through tet
        bool cutTetCell
        (
            const scalar target,
            const label faceOrig,
            const bool ow
        );

        //- Calculate local addressing for multi-cut
        inline void calcAddressing
        (
            const MPLICcellStorage& cellInfo
        );

        //- Append face area vectors and centers to cache
        inline void appendSfCf
        (
            const vector& Sf,
            const vector& Cf,
            const scalar magSf,
            const bool own = true
        );

        //- Calculating surface vector of unordered edges
        inline bool cutStatusCalcSf();

        //- Calculate cut surface area vector
        inline vector calcCutSf() const;

        //- Calculating surface center of unordered edges
        inline vector calcCutCf(const vector& cutSf) const;

        //- Calculate phiU
        inline void phiU
        (
            const pointField& points,
            const faceList& faces,
            const labelList& cFaces,
            const vectorField& pointsU
        );

        //- Resize and set all the cell face fields to 0
        inline void resetFaceFields(const label nFaces);

        //- Compute surface interpolation from area ratio
        inline void calcAlphaf(const UIndirectList<scalar>& magSfs);

        //- Compute surface interpolation from flux ratio
        inline void calcAlphaUf();

        //- Clear storage
        inline void clear();

        //- Clear single cut storage
        inline void clearOneCut();


public:

    // Constructors

        //- Construct for given interpolation and PLIC type
        MPLICcell
        (
            const bool unweighted = true,
            const bool multiCut = true
        );


    // Member functions

        //- Match cell volume ratios
        bool matchAlpha
        (
            const MPLICcellStorage& cellInfo
        );

        //- Return face volume fraction
        inline const DynamicList<scalar>& alphaf() const;

        //- Return cut normal
        inline const vector& cutNormal() const;

        //- Return volume fraction corresponding to the cut
        inline scalar cutAlpha() const;

        //- Return sub-cell volume
        inline scalar subCellVolume() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MPLICcellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
