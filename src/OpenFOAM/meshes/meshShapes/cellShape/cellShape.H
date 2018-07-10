/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::cellShape

Description
    An analytical geometric cellShape.

    The optional collapse functionality changes the cellModel to the
    correct type after removing any duplicate points.

SourceFiles
    cellShapeI.H
    cellShape.C
    cellShapeIO.C
    cellShapeEqual.C

\*---------------------------------------------------------------------------*/

#ifndef cellShape_H
#define cellShape_H

#include "pointField.H"
#include "labelList.H"
#include "cellModel.H"
#include "autoPtr.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cell;

// Forward declaration of friend functions and operators

class cellShape;
bool operator==(const cellShape&, const cellShape&);
Istream& operator>>(Istream&, cellShape&);
Ostream& operator<<(Ostream&, const cellShape&);


/*---------------------------------------------------------------------------*\
                          Class cellShape Declaration
\*---------------------------------------------------------------------------*/

class cellShape
:
    public labelList
{
    // Private data

        //- Access to the cellShape's model
        const cellModel *m;


public:

    // Constructors

        //- Construct null
        inline cellShape();

        //- Construct from components
        inline cellShape
        (
            const cellModel&,
            const labelList&,
            const bool doCollapse = false
        );

        //- Construct from components
        inline cellShape
        (
            const word& model,
            const labelList&,
            const bool doCollapse = false
        );

        //- Construct from Istream
        inline cellShape(Istream& is);

        //- Clone
        inline autoPtr<cellShape> clone() const;


    // Member Functions

        //- Return the points corresponding to this cellShape
        inline pointField points(const pointField& meshPoints) const;

        //- Model reference
        inline const cellModel& model() const;

        //- Mesh face labels of this cell (in order of model)
        inline labelList meshFaces(const faceList& allFaces, const cell&) const;

        //- Mesh edge labels of this cell (in order of model)
        inline labelList meshEdges
        (
            const edgeList& allEdges,
            const labelList&
        ) const;

        //- Faces of this cell
        inline faceList faces() const;

        //- Collapsed faces of this cell
        inline faceList collapsedFaces() const;

        //- Number of faces
        inline label nFaces() const;

        //- Edges of this cellShape
        inline edgeList edges() const;

        //- Number of edges
        inline label nEdges() const;

        //- Number of points
        inline label nPoints() const;

        //- Centroid of the cell
        inline point centre(const pointField&) const;

        //- Return info proxy.
        //  Used to print token information to a stream
        Foam::InfoProxy<cellShape> info() const
        {
            return *this;
        }

        //- Scalar magnitude
        inline scalar mag(const pointField&) const;

        //- Collapse shape to correct one after removing duplicate vertices
        void collapse();


    // Friend Operators

        friend bool operator==(const cellShape&, const cellShape&);


    // IOstream operators

        friend Istream& operator>>(Istream&, cellShape&);
        friend Ostream& operator<<(Ostream&, const cellShape&);
};


template<>
Ostream& operator<<(Ostream& os, const InfoProxy<cellShape>& ip);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellShapeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
