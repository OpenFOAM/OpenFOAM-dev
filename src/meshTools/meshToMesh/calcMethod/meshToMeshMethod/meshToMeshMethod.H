/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::meshToMeshMethod

Description
    Base class for mesh-to-mesh calculation methods

SourceFiles
    meshToMeshMethod.C

\*---------------------------------------------------------------------------*/

#ifndef meshToMeshMethod_H
#define meshToMeshMethod_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class meshToMeshMethod Declaration
\*---------------------------------------------------------------------------*/

class meshToMeshMethod
{
    // Private static data

        //- Map from mesh-to-mesh method to patch-to-patch method
        static HashTable<word> meshToMeshToPatchToPatchMethod_;


protected:

    // Protected data

        //- Reference to the source mesh
        const polyMesh& src_;

        //- Reference to the target mesh
        const polyMesh& tgt_;


    // Protected Member Functions

        //- Return src cell IDs for the overlap region
        labelList maskCells() const;

        //- Append target cell neighbour cells to cellIDs list
        void appendNbrCells
        (
            const label tgtCelli,
            const polyMesh& mesh,
            const DynamicList<label>& visitedTgtCells,
            DynamicList<label>& nbrTgtCellIDs
        ) const;

        //- Initialise the arrays used to store the intersection result. Return
        //  whether or not an intersection is possible.
        bool initialise
        (
            labelListList& srcToTgtAddr,
            scalarListList& srcToTgtWght,
            labelListList& tgtToTgtAddr,
            scalarListList& tgtToTgtWght
        ) const;

        //- Write the connectivity for debugging
        void writeConnectivity
        (
            const polyMesh& mesh1,
            const polyMesh& mesh2,
            const labelListList& mesh1ToMesh2Addr
        ) const;


public:

    //- Run-time type information
    TypeName("meshToMeshMethod");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        meshToMeshMethod,
        components,
        (
            const polyMesh& src,
            const polyMesh& tgt
        ),
        (src, tgt)
    );


    // Constructors

        //- Construct from source and target meshes
        meshToMeshMethod(const polyMesh& src, const polyMesh& tgt);

        //- Disallow default bitwise copy construction
        meshToMeshMethod(const meshToMeshMethod&) = delete;


    //- Selector
    static autoPtr<meshToMeshMethod> New
    (
        const word& methodName,
        const polyMesh& src,
        const polyMesh& tgt
    );


    //- Destructor
    virtual ~meshToMeshMethod();


    // Member Functions

        //- Get the corresponding patch-to-patch method
        virtual const word& patchToPatchMethod() const = 0;


        // Evaluate

            //- Calculate addressing and weights
            virtual scalar calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                labelListList& tgtToTgtAddr,
                scalarListList& tgtToTgtWght
            ) = 0;

            //- Normalise the weights for a given mesh
            virtual void normalise
            (
                const polyMesh& srcMesh,
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght
            ) const = 0;


        // Access

            //- Return const access to the source mesh
            inline const polyMesh& src() const
            {
                return src_;
            }

            //- Return const access to the target mesh
            inline const polyMesh& tgt() const
            {
                return tgt_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMeshMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
