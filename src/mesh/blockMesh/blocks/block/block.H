/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    // Private Data

        //- List of points
        pointField points_;

        //- Boundary patches
        FixedList<List<FixedList<label, 4>>, 6> boundaryPatches_;


    // Private Member Functions

        //- Creates vertices for cells filling the block
        void createPoints();

        //- Creates boundary patch faces for the block
        void createBoundary();


public:

    //- Runtime type information
    TypeName("block");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            block,
            Istream,
            (
                const dictionary& dict,
                const label index,
                const pointField& vertices,
                const blockEdgeList& edges,
                const blockFaceList& faces,
                Istream& is
            ),
            (dict, index, vertices, edges, faces, is)
        );


    // Constructors

        //- Construct from components with Istream
        block
        (
            const dictionary& dict,
            const label index,
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream& is
        );

        //- Construct from a block definition
        block(const blockDescriptor&);

        //- Disallow default bitwise copy construction
        block(const block&) = delete;

        //- Clone
        autoPtr<block> clone() const
        {
            NotImplemented;
            return autoPtr<block>(nullptr);
        }

        //- New function which constructs and returns pointer to a block
        static autoPtr<block> New
        (
            const dictionary& dict,
            const label index,
            const pointField& points,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream&
        );

        //- Class used for the read-construction of
        //  PtrLists of blocks
        class iNew
        {
            const dictionary& dict_;
            const pointField& points_;
            const blockEdgeList& edges_;
            const blockFaceList& faces_;
            mutable label index_;

        public:

            iNew
            (
                const dictionary& dict,
                const pointField& points,
                const blockEdgeList& edges,
                const blockFaceList& faces
            )
            :
                dict_(dict),
                points_(points),
                edges_(edges),
                faces_(faces),
                index_(0)
            {}

            autoPtr<block> operator()(Istream& is) const
            {
                return block::New(dict_, index_++, points_, edges_, faces_, is);
            }
        };


    //- Destructor
    virtual ~block()
    {}


    // Member Functions

        // Access

            //- Return the points for filling the block
            inline const pointField& points() const;

            //- Return the cells for filling the block
            List<FixedList<label, 8>> cells() const;

            //- Return the boundary patch faces for the block
            inline const FixedList<List<FixedList<label, 4>>, 6>&
                boundaryPatches() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const block&) = delete;


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
