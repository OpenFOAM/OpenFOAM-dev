/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::pointFieldDecomposer

Description
    Point field decomposer.

SourceFiles
    pointFieldDecomposer.C
    pointFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef pointFieldDecomposer_H
#define pointFieldDecomposer_H

#include "pointMesh.H"
#include "pointPatchFieldMapperPatchRef.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pointFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class pointFieldDecomposer
{

public:

        //- Point patch field decomposer class
        class patchFieldDecomposer
        :
            public pointPatchFieldMapperPatchRef
        {
            // Private Data

                labelList directAddressing_;

                //- Does map contain any unmapped values
                bool hasUnmapped_;

        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const pointPatch& completeMeshPatch,
                    const pointPatch& procMeshPatch,
                    const labelList& directAddr
                );


            // Member Functions

                label size() const
                {
                    return directAddressing_.size();
                }

                bool direct() const
                {
                    return true;
                }

                bool hasUnmapped() const
                {
                    return hasUnmapped_;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


private:

    // Private Data

        //- Reference to complete mesh
        const pointMesh& completeMesh_;

        //- Reference to processor mesh
        const pointMesh& procMesh_;

        //- Reference to point addressing
        const labelList& pointAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        List<patchFieldDecomposer*> patchFieldDecomposerPtrs_;


public:

    // Constructors

        //- Construct from components
        pointFieldDecomposer
        (
            const pointMesh& completeMesh,
            const pointMesh& procMesh,
            const labelList& pointAddressing,
            const labelList& boundaryAddressing
        );

        //- Disallow default bitwise copy construction
        pointFieldDecomposer(const pointFieldDecomposer&) = delete;


    //- Destructor
    ~pointFieldDecomposer();


    // Member Functions

        //- Decompose point field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        decomposeField
        (
            const GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
