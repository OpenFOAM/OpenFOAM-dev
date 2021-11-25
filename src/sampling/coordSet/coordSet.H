/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::coordSet

Description
    Holds list of sampling positions

SourceFiles
    coordSet.C

\*---------------------------------------------------------------------------*/

#ifndef coordSet_H
#define coordSet_H

#include "pointField.H"
#include "word.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class coordSet Declaration
\*---------------------------------------------------------------------------*/

class coordSet
{
public:

    // Public data types

        //- Enumeration defining the output format for coordinates
        enum class axisType
        {
            XYZ,
            X,
            Y,
            Z,
            DISTANCE,
            DEFAULT
        };


        //- String representation of axis enums
        static const NamedEnum<axisType, 6> axisTypeNames_;


protected:

    // Protected Data

        //- Connected segments
        labelList segments_;

        //- Name of the positions
        word positionName_;

        //- Point positions
        autoPtr<pointField> positions_;

        //- Name of the distances
        word distanceName_;

        //- Scalar distances
        autoPtr<scalarField> distances_;

        //- Axis
        axisType axis_;


public:

    // Constructors

        //- Construct null
        coordSet();

        //- Construct from components
        coordSet
        (
            const labelList& segments,
            const word& positionName = word::null,
            const pointField& positions = pointField::null(),
            const word& distanceName = word::null,
            const scalarField& distances = scalarField::null(),
            const word& axis = axisTypeNames_[axisType::DEFAULT]
        );

        //- Construct from positions
        coordSet
        (
            const bool contiguous,
            const word& positionName,
            const pointField& positions,
            const word& axis = axisTypeNames_[axisType::DEFAULT]
        );

        //- Construct from distances
        coordSet
        (
            const bool contiguous,
            const word& distanceName,
            const scalarField& distances,
            const word& axis = axisTypeNames_[axisType::DEFAULT]
        );


    // Member Functions

        //- Return the size
        inline label size() const
        {
            return segments_.size();
        }

        //- Return the segments
        inline const labelList& segments() const
        {
            return segments_;
        }

        //- Return the axis name
        inline word axis() const
        {
            return axisTypeNames_[axis_];
        }

        //- Is the coordinate axis a scalar?
        bool hasScalarAxis() const;

        //- Is the coordinate axis a point?
        bool hasPointAxis() const;

        //- Get scalar coordinate (axis is x, y, z or distance)
        scalar scalarCoord(const label index) const;

        //- Get scalar coordinates (axis is x, y, z or distance)
        tmp<scalarField> scalarCoords() const;

        //- Return the name of the scalar coordinates
        word scalarName() const;

        //- Get vector coordinate (axis is xyz)
        point pointCoord(const label index) const;

        //- Get vector coordinate (axis is xyz)
        tmp<pointField> pointCoords() const;

        //- Return the name of the point coordinates
        word pointName() const;

        //- Return a list of isolated vertices. These are the points that are
        //  not adjacent to any points in the same segment.
        labelList vertices() const;

        //- Return a list of edges. These are adjacent pairs of points which
        //  are in the same segment.
        labelPairList edges() const;

        //- Return a list of lines. These are lists of points which are in the
        //  same segment.
        labelListList lines() const;

        //- Combine coordinate sets onto the master. Return both the combined
        //  coordinate set, and an ordering to be used for gathering associated
        //  fields
        Tuple2<coordSet, labelList> gather() const;

        //- Combine a field using the ordering obtained from the coordinate set
        //  gather operation
        template<class Type>
        static tmp<Field<Type>> gather
        (
            const Field<Type>& values,
            const labelList& order
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coordSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
