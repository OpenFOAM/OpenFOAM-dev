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
    Foam::blockDescriptor

Description
    Takes the description of the block and the list of curved edges and
    creates a list of points on edges together with the weighting factors

SourceFiles
    blockDescriptor.C
    blockDescriptorEdges.C

\*---------------------------------------------------------------------------*/

#ifndef blockDescriptor_H
#define blockDescriptor_H

#include "cellShape.H"
#include "pointField.H"
#include "scalarList.H"
#include "curvedEdgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class blockMesh;
class blockDescriptor;
Ostream& operator<<(Ostream&, const blockDescriptor&);

/*---------------------------------------------------------------------------*\
                       Class blockDescriptor Declaration
\*---------------------------------------------------------------------------*/

class blockDescriptor
{
    // Private data

        //- Reference to point field defining the block mesh
        const pointField& blockPointField_;

        //- Reference to a list of curved edges
        const curvedEdgeList& curvedEdges_;

        //- Block shape
        cellShape blockShape_;

        //- The number of cells in the i,j,k directions
        Vector<label> meshDensity_;

        //- Block edge points
        List< List<point> > edgePoints_;

        //- Block edge weighting factors
        scalarListList edgeWeights_;

        //- Expansion ratios in all directions
        scalarList expand_;

        //- Name of the zone (empty string if none)
        word zoneName_;

    // Private Member Functions

        //- Set the points/weights for all edges
        void makeBlockEdges();

        //- Set the edge points/weights
        void setEdge(label edgeI, label start, label end, label dim);


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const blockDescriptor&);


public:

    // Constructors

        //- Construct from components. Optional cellSet/zone name.
        blockDescriptor
        (
            const cellShape&,
            const pointField& blockPointField,
            const curvedEdgeList&,
            const Vector<label>& meshDensity,
            const UList<scalar>& expand,
            const word& zoneName = ""
        );

        //- Construct from Istream
        blockDescriptor
        (
            const pointField& blockPointField,
            const curvedEdgeList&,
            Istream&
        );

        //- Clone
        autoPtr<blockDescriptor> clone() const
        {
            notImplemented("blockDescriptor::clone()");
            return autoPtr<blockDescriptor>(NULL);
        }


    //- Destructor
    ~blockDescriptor();


    // Member Functions

        // Access

        //- Reference to point field defining the block mesh
        const pointField& blockPointField() const;

        //- Return the block shape
        const cellShape& blockShape() const;

        //- Return the block points along each edge
        const List< List<point> >& blockEdgePoints() const;

        //- Return the weightings along each edge
        const scalarListList& blockEdgeWeights() const;

        //- Return the mesh density (number of cells) in the i,j,k directions
        const Vector<label>& meshDensity() const;

        //- Return the (optional) zone name
        const word& zoneName() const;


        //- Return the number of points
        label nPoints() const;

        //- Return the number of cells
        label nCells() const;

        //- Return block point at local label i
        const point& blockPoint(const label i) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const blockDescriptor&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
