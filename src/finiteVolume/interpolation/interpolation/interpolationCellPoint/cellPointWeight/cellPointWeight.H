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
    Foam::cellPointWeight

Description
    Foam::cellPointWeight

SourceFiles
    cellPointWeight.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointWeight_H
#define cellPointWeight_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class cellPointWeight Declaration
\*---------------------------------------------------------------------------*/

class cellPointWeight
{
protected:

    // Protected data

       //- Cell index
       const label cellI_;

       //- Weights applied to tet vertices
       List<scalar> weights_;

       //- Face vertex indices
       List<label> faceVertices_;


    // Protected Member Functions

        void findTetrahedron
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI
        );

        void findTriangle
        (
            const polyMesh& mesh,
            const vector& position,
            const label faceI
        );


public:

    //- Debug switch
    static int debug;

    //- Tolerance used in calculating barycentric co-ordinates
    //  (applied to normalised values)
    static scalar tol;


    // Constructors

        //- Construct from components
        cellPointWeight
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label faceI = -1
        );


    // Member Functions

        //- Cell index
        inline label cell() const
        {
            return cellI_;
        }

        //- interpolation weights
        inline const List<scalar>& weights() const
        {
            return weights_;
        }

        //- interpolation addressing for points on face
        inline const List<label>& faceVertices() const
        {
            return faceVertices_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
