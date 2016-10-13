/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::block

Description
    Creates a single block of cells from point coordinates, numbers of
    cells in each direction and an expansion ratio.

Note
    The vertices and cells for filling the block are demand-driven.

SourceFiles
    block.C
    blockCreate.C

\*---------------------------------------------------------------------------*/

#ifndef block_H
#define block_H

#include "pointField.H"
#include "labelList.H"

#include "blockDescriptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class block;
Ostream& operator<<(Ostream&, const block&);

/*---------------------------------------------------------------------------*\
                           Class block Declaration
\*---------------------------------------------------------------------------*/

class block
:
    public blockDescriptor
{
    // Private data

        //- List of points
        pointField points_;

        //- Boundary patches
        FixedList<List<FixedList<label, 4>>, 6> boundaryPatches_;


    // Private Member Functions

        //- Creates vertices for cells filling the block
        void createPoints();

        //- Creates boundary patch faces for the block
        void createBoundary();

        //- Disallow default bitwise copy construct
        block(const block&);

        //- Disallow default bitwise assignment
        void operator=(const block&);


public:

    // Constructors

        //- Construct from components with Istream
        block
        (
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream& is
        );

        //- Construct from a block definition
        block(const blockDescriptor&);

        //- Clone
        autoPtr<block> clone() const
        {
            NotImplemented;
            return autoPtr<block>(nullptr);
        }

        //- Class used for the read-construction of
        //  PtrLists of blocks
        class iNew
        {
            const pointField& points_;
            const blockEdgeList& edges_;
            const blockFaceList& faces_;

        public:

            iNew
            (
                const pointField& points,
                const blockEdgeList& edges,
                const blockFaceList& faces
            )
            :
                points_(points),
                edges_(edges),
                faces_(faces)
            {}

            autoPtr<block> operator()(Istream& is) const
            {
                return autoPtr<block>(new block(points_, edges_, faces_, is));
            }
        };


    // Member Functions

        // Access

            //- Return the points for filling the block
            inline const pointField& points() const;

            //- Return the cells for filling the block
            List<FixedList<label, 8>> cells() const;

            //- Return the boundary patch faces for the block
            inline const FixedList<List<FixedList<label, 4>>, 6>&
                boundaryPatches() const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const block&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
