/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "pointPatchFieldMapperPatchRef.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pointFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class pointFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        const pointMesh& mesh_;

        //- List of processor meshes
        const PtrList<pointMesh>& procMeshes_;

        //- List of processor point addressing lists
        const PtrList<labelIOList>& pointProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;

        //- Point patch addressing
        labelListListList patchPointAddressing_;

        //- Number of fields reconstructed
        label nReconstructed_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointFieldReconstructor(const pointFieldReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const pointFieldReconstructor&);


public:

        class pointPatchFieldReconstructor
        :
            public pointPatchFieldMapper
        {
            label size_;

        public:

            // Constructors

                //- Construct given size
                pointPatchFieldReconstructor(const label size)
                :
                    size_(size)
                {}


            // Member functions

                label size() const
                {
                    return size_;
                }

                bool direct() const
                {
                    return true;
                }

                bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return labelUList::null();
                }
        };


    // Constructors

        //- Construct from components
        pointFieldReconstructor
        (
            const pointMesh& mesh,
            const PtrList<pointMesh>& procMeshes,
            const PtrList<labelIOList>& pointProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Return number of fields reconstructed
        label nReconstructed() const
        {
            return nReconstructed_;
        }

        //- Reconstruct field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> >
        reconstructField(const IOobject& fieldIoObject);

        //- Reconstruct and write all fields
        template<class Type>
        void reconstructFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
