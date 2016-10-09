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
    Foam::blockDescriptor

Description
    Takes the description of the block and the list of curved edges and
    creates a list of points on edges together with the weighting factors

    For a given block, the correspondence between the ordering of vertex labels
    and face labels is shown below.  For vertex numbering in the sequence 0 to 7
    (block, centre): faces 0 (f0) and 1 are left and right, respectively; faces
    2 and 3 are bottom and top; and faces 4 and 5 are front the back:
    \verbatim
             4 ---- 5
        f3   |\     |\   f5
         |   | 7 ---- 6    \
         |   0 |--- 1 |     \
         |    \|     \|      f4
        f2     3 ---- 2

             f0 ----- f1
     \endverbatim

SourceFiles
    blockDescriptor.C
    blockDescriptorEdges.C

\*---------------------------------------------------------------------------*/

#ifndef blockDescriptor_H
#define blockDescriptor_H

#include "cellShape.H"
#include "pointField.H"
#include "scalarList.H"
#include "blockEdgeList.H"
#include "gradingDescriptors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class blockDescriptor;
Ostream& operator<<(Ostream&, const blockDescriptor&);

/*---------------------------------------------------------------------------*\
                       Class blockDescriptor Declaration
\*---------------------------------------------------------------------------*/

class blockDescriptor
{
    // Private data

        //- Reference to point field defining the block mesh
        const pointField& vertices_;

        //- Reference to a list of block edges
        const blockEdgeList& edges_;

        //- Block shape
        cellShape blockShape_;

        //- The number of cells in the i,j,k directions
        Vector<label> density_;

        //- Expansion ratios in all directions
        List<gradingDescriptors> expand_;

        //- Name of the zone (empty word if none)
        word zoneName_;


    // Private Member Functions

        //- Check block has outward-pointing faces
        void check(const Istream& is);

        //- Calculate the points and weights for the specified edge
        void edgePointsWeights
        (
            List<point> (&edgePoints)[12],
            scalarList (&edgeWeights)[12],
            const label edgei,
            const label start,
            const label end,
            const label dim
        ) const;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const blockDescriptor&);


public:

    // Constructors

        //- Construct from components. Optional cellSet/zone name.
        blockDescriptor
        (
            const cellShape&,
            const pointField& vertices,
            const blockEdgeList&,
            const Vector<label>& density,
            const UList<gradingDescriptors>& expand,
            const word& zoneName = ""
        );

        //- Construct from Istream
        blockDescriptor
        (
            const pointField& vertices,
            const blockEdgeList&,
            Istream&
        );


    //- Destructor
    ~blockDescriptor();


    // Member Functions

        // Access

        //- Return the block shape
        inline const cellShape& blockShape() const;

        //- Return the mesh density (number of cells) in the i,j,k directions
        inline const Vector<label>& density() const;

        //- Return the (optional) zone name
        inline const word& zoneName() const;

        //- Return the number of points
        inline label nPoints() const;

        //- Return the number of cells
        inline label nCells() const;

        //- Return block point for local label i
        inline const point& blockPoint(const label i) const;

        //- Calculate the points and weights for all edges
        void edgesPointsWeights
        (
            List<point> (&edgePoints)[12],
            scalarList (&edgeWeights)[12]
        ) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const blockDescriptor&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockDescriptorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
