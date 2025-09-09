/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::extrudeModels::path

Description
    Extrudes by transforming points along an curve defined as an edgeMesh

    Uses rotation minimising frames (RMF) algorithm, which is based on
    the Frenet-Serret equations, and described in:
    \verbatim
        Wang, W., Jüttler, B., Zheng, D., & Liu, Y. (2008).
        Computation of rotation minimizing frames.
        ACM Transactions on Graphics (TOG), 27(1), 1-18.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef path_extrudeModel_H
#define path_extrudeModel_H

#include "extrudeModel.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class edgeMesh;

namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                      Class path Declaration
\*---------------------------------------------------------------------------*/

class path
:
    public extrudeModel
{
    // Private Data

        //- Coordinates defining the extrusion curve
        autoPtr<edgeMesh> eMeshPtr_;

        //- Cumulative point distances along the extrusion curve
        scalarField distances_;

        //- Extrude direction for each frame along the extrusion curve
        vectorField directions_;

        //- Normals to each direction frame
        vectorField normals_;

        //- Transposed orthonormal rotation tensor at the first layer
        tensor R0T_;


    // Private Member Functions

        //- Find cumulative point span which contains the given distance
        label findFrameIndex(const scalar distance) const;


    // Private Static Member Functions

        //- Cartesian axis perpendicular to the supplied direction
        static vector orthogonal(const vector& dirn);

        //- Frenet-Serret rotational tensor
        static tensor orthonormalBasis(const vector& T, const vector& N);


public:

    //- Runtime type information
    TypeName("path");


    // Constructors

        //- Construct from dictionary
        path(const dictionary& dict);


    //- Destructor
    virtual ~path();


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
