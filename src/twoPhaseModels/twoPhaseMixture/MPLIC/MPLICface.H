/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::MPLICface

Description
    Class that deals with cutting faces based on face point values and target
    value.

SourceFiles
    MPLICface.C

\*---------------------------------------------------------------------------*/

#ifndef MPLICface_H
#define MPLICface_H

#include "face.H"
#include "pointField.H"
#include "DynamicList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MPLICface Declaration
\*---------------------------------------------------------------------------*/

class MPLICface
{
    // Private Member variables

        //- Store all the cut points
        DynamicList<point> cutPoints_;

        //- Store points of sub face
        DynamicList<point> subPoints_;

        //- Labels of cut edges
        DynamicList<label> cutEdges_;

        //- Store velocity point values
        DynamicList<vector> subPointsU_;

        //- Keep truck of face orientation
        bool flipped_;

        //- Select unweighted interpolation if true
        // velocity flux corrected if false
        const bool unweighted_;


    public:

        // Constructors

            //- Construct empty
            MPLICface(const bool unweighted = true);

            // - Destructor
            ~MPLICface()
            {}


        // Member functions

            //- Function to cut for multi cut
            //  Returns:
            //  - 0: no cut
            //  - 1: cut
            label cutFace
            (
                const labelList& f,
                const labelList& faceEdges,
                const pointField& points,
                const boolList& isEdgeCutOld,
                boolList& isEdgeCut,
                label& faceEdgei,
                const UList<scalar>& pointsAlpha,
                const UList<vector>& pointsU,
                const label facei,
                const scalar target,
                const bool ow
            );

            //- Cut the face and return the type of cut
            //  Returns:
            //  - -1: multi cut
            //  -  0: no cut
            //  - +1: single cut
            label cutFace
            (
                const UList<label>& f,
                const UList<point>& points,
                const UList<scalar>& pointsAlpha,
                const UList<vector>& pointsU,
                const scalar target,
                const bool ow
            );

            //- Calculate and return alphaPhiU
            inline scalar alphaPhiU() const;

            //- Calculate and return alphaPhiU
            template<class VectorList, class PointList>
            inline scalar alphaPhiU
            (
                const VectorList& pointsU,
                const PointList& points
            ) const;

            //- Calculate and return alphaPhiU
            template<class VectorList, class PointList>
            inline scalar alphaPhiU
            (
                const VectorList& pointsU,
                const PointList& points,
                const labelList& f
            ) const;

            //- Access to cut points
            inline const DynamicList<point>& cutPoints() const;

            //- Access to submerged face points
            inline const DynamicList<point>& subPoints() const;

            //- Access to cut edges
            inline const DynamicList<label>& cutEdges() const;

            //- Return subface surface area vector
            inline const vector Sf() const;

            //- Return subface centre
            inline const vector Cf(const vector& area) const;

            //- Return interpolated U values
            inline const DynamicList<vector>& U() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MPLICfaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
