/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::leastSquaresVectors

Description
    Least-squares gradient scheme vectors

SourceFiles
    leastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresVectors_H
#define leastSquaresVectors_H

#include "DemandDrivenMeshObject.H"
#include "fvMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresVectors
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        leastSquaresVectors
    >
{
    // Private Data

        //- Least-squares gradient vectors
        surfaceVectorField pVectors_;
        surfaceVectorField nVectors_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void calcLeastSquaresVectors();


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        leastSquaresVectors
    >;

    // Protected Constructors

        //- Construct given an fvMesh
        explicit leastSquaresVectors(const fvMesh&);


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresVectors");


    //- Destructor
    virtual ~leastSquaresVectors();


    // Member Functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const
        {
            return pVectors_;
        }

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const
        {
            return nVectors_;
        }

        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
