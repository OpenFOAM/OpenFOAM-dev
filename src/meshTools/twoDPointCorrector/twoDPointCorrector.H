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
    Foam::twoDPointCorrector

Description
    Class applies a two-dimensional correction to mesh motion point field.

    The correction guarantees that the mesh does not get twisted during motion
    and thus introduce a third dimension into a 2-D problem.

    The operation is performed by looping through all edges approximately
    normal to the plane and enforcing their orthogonality onto the plane by
    adjusting points on their ends.

SourceFiles
    twoDPointCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef twoDPointCorrector_H
#define twoDPointCorrector_H

#include "MeshObject.H"
#include "pointField.H"
#include "labelList.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class twoDPointCorrector Declaration
\*---------------------------------------------------------------------------*/

class twoDPointCorrector
:
    public MeshObject<polyMesh, UpdateableMeshObject, twoDPointCorrector>
{
    // Private Data

        //- Is 2D correction required, i.e. is the mesh
        bool required_;

        //- 2-D plane unit normal
        mutable vector* planeNormalPtr_;

        //- Indices of edges normal to plane
        mutable labelList* normalEdgeIndicesPtr_;

        //- Flag to indicate a wedge geometry
        mutable bool isWedge_;

        //- Wedge axis (if wedge geometry)
        mutable vector wedgeAxis_;

        //- Wedge angle (if wedge geometry)
        mutable scalar wedgeAngle_;


    // Private Member Functions

        //- Calculate addressing
        void calcAddressing() const;

        //- Clear addressing
        void clearAddressing() const;

        //- Snap a point to the wedge patch(es)
        void snapToWedge(const vector& n, const point& A, point& p) const;


    // Static Data Members

        //- Edge orthogonality tolerance
        static const scalar edgeOrthogonalityTol;


public:

    // Declare name of the class and its debug switch
    ClassName("twoDPointCorrector");


    // Constructors

        //- Construct from components
        twoDPointCorrector(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        twoDPointCorrector(const twoDPointCorrector&) = delete;


    //- Destructor
    ~twoDPointCorrector();


    // Member Functions

        //- Is 2D correction required, i.e. is the mesh a wedge or slab
        bool required() const
        {
            return required_;
        }

        //- Return plane normal
        const vector& planeNormal() const;

        //- Return indices of normal edges.
        const labelList& normalEdgeIndices() const;

        //- Return direction normal to plane
        direction normalDir() const;

        //- Correct motion points
        void correctPoints(pointField& p) const;

        //- Correct motion displacements
        void correctDisplacement(const pointField& p, vectorField& disp) const;

        //- Update topology
        void updateMesh(const mapPolyMesh&);

        //- Correct weighting factors for moving mesh.
        bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const twoDPointCorrector&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
