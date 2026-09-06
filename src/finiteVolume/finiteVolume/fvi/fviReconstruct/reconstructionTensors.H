/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::reconstructionTensors

Description
    Flux reconstruction tensors

SourceFiles
    reconstructionTensors.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructionTensors_H
#define reconstructionTensors_H

#include "DemandDrivenMeshObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvi
{

/*---------------------------------------------------------------------------*\
                   Class reconstructionTensors Declaration
\*---------------------------------------------------------------------------*/

class reconstructionTensors
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        reconstructionTensors
    >
{
    // Private Data

        //- Reconstruction tensors
        VolInternalField<symmTensor> tensors_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void calcReconstructTensors();


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        reconstructionTensors
    >;

    // Protected Constructors

        //- Construct given an fvMesh
        explicit reconstructionTensors(const fvMesh&);


public:

    // Declare name of the class and its debug switch
    TypeName("reconstructionTensors");


    //- Destructor
    virtual ~reconstructionTensors();


    // Member Functions

        //- Return reference to the reconstruction tensors
        const VolInternalField<symmTensor>& tensors() const
        {
            return tensors_;
        }

        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
