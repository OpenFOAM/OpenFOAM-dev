/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::pointFieldReconstructor

Description
    Point field reconstructor.

SourceFiles
    pointFieldReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef pointFieldReconstructor_H
#define pointFieldReconstructor_H

#include "pointMesh.H"
#include "pointFields.H"
#include "pointPatchFieldMapper.H"
#include "setSizePointPatchFieldMapper.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class pointFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class pointFieldReconstructor
{
    // Private Data

        //- Reconstructed mesh reference
        const pointMesh& completeMesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- List of processor point addressing lists
        const labelListList& pointProcAddressing_;

        //- Point patch addressing
        labelListListList patchPointAddressing_;

        //- Number of fields reconstructed
        label nReconstructed_;


public:

    // Constructors

        //- Construct from components
        pointFieldReconstructor
        (
            const pointMesh& mesh,
            const PtrList<fvMesh>& procMeshes,
            const labelListList& pointProcAddressing
        );

        //- Disallow default bitwise copy construction
        pointFieldReconstructor(const pointFieldReconstructor&) = delete;


    // Member Functions

        //- Return number of fields reconstructed
        label nReconstructed() const
        {
            return nReconstructed_;
        }

        //- Reconstruct field
        template<class Type>
        tmp<PointField<Type>>
        reconstructField(const IOobject& fieldIoObject);

        //- Reconstruct and write all fields
        template<class Type>
        void reconstructFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointFieldReconstructor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
