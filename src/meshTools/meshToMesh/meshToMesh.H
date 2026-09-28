/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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
    Foam::meshToMesh

Description
    Class to calculate interpolative addressing and weights between the cells
    and patches of two overlapping meshes

SourceFiles
    meshToMesh.C

\*---------------------------------------------------------------------------*/

#ifndef meshToMesh_H
#define meshToMesh_H

#include "polyMesh.H"
#include "cellsToCells.H"
#include "cellsToCellsExtrapolation.H"
#include "patchToPatch.H"
#include "patchToPatchExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshToMesh Declaration
\*---------------------------------------------------------------------------*/

class meshToMesh
{
    // Private Data

        //- Reference to the source mesh
        const polyMesh& srcMesh_;

        //- Reference to the target mesh
        const polyMesh& tgtMesh_;

        //- Interpolation engine between source and target cells
        autoPtr<cellsToCells> cellsInterpolation_;

        //- Cells extrapolation engine type
        const word cellsExtrapolationType_;

        //- Extrapolation engine for the source cells
        mutable autoPtr<cellsToCellsExtrapolation> srcCellsExtrapolation_;

        //- Extrapolation engine for the target cells
        mutable autoPtr<cellsToCellsExtrapolation> tgtCellsExtrapolation_;

        //- List of corresponding source and target patch indices
        List<labelPair> patchIndices_;

        //- List of interpolation engines between source and target patches
        PtrList<patchToPatch> patchInterpolations_;

        //- Patch extrapolation engine type
        const word patchExtrapolationType_;

        //- List of extrapolation engines for the source patches
        mutable PtrList<patchToPatchExtrapolation> srcPatchExtrapolations_;

        //- List of extrapolation engines for the source patches
        mutable PtrList<patchToPatchExtrapolation> tgtPatchExtrapolations_;


protected:

    // Protected Member Functions

        // Access

            //- Return the interpolation engine between source and target cells
            inline const cellsToCells& cellsInterpolation() const;

            //- Return the extrapolation engine for the source cells
            inline const cellsToCellsExtrapolation&
                srcCellsExtrapolation() const;

            //- Return the extrapolation engine for the target cells
            inline const cellsToCellsExtrapolation&
                tgtCellsExtrapolation() const;

            //- Return the list of corresponding source and target patch indices
            inline const List<labelPair>& patchIndices() const;

            //- Return the interpolation engine between a source and a target
            //  patch
            inline const patchToPatch& patchInterpolation
            (
                const label i
            ) const;

            //- Return the extrapolation engine for a source patch
            inline const patchToPatchExtrapolation& srcPatchExtrapolation
            (
                const label i
            ) const;

            //- Return the extrapolation engine for a target patch
            inline const patchToPatchExtrapolation& tgtPatchExtrapolation
            (
                const label i
            ) const;


public:

    //- Run-time type information
    TypeName("meshToMesh");


    // Constructors

        //- Construct from source and target meshes. If a patchMap is supplied,
        //  then interpolate between the specified patches. If not, then assume
        //  a consistent mesh with consistently named patches and interpolate
        //  1-to-1 between patches with the same name.
        meshToMesh
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            const word& cellsInterpolationType,
            const word& cellsExtrapolationType,
            const word& patchInterpolationType,
            const word& patchExtrapolationType,
            const HashTable<word>& patchMap = NullObjectRef<HashTable<word>>()
        );

        //- As above but with default extrapolation methods
        meshToMesh
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            const word& cellsInterpolationType,
            const word& patchInterpolationType,
            const HashTable<word>& patchMap = NullObjectRef<HashTable<word>>()
        );

        //- As above but with the same interpolation type for both cells and
        //  patches and default extrapolation methods
        meshToMesh
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            const word& interpolationType,
            const HashTable<word>& patchMap = NullObjectRef<HashTable<word>>()
        );

        //- Disallow default bitwise copy construction
        meshToMesh(const meshToMesh&) = delete;


    //- Destructor
    virtual ~meshToMesh();


    // Member Functions

        // Access

            //- Return const access to the source mesh
            inline const polyMesh& srcMesh() const;

            //- Return const access to the target mesh
            inline const polyMesh& tgtMesh() const;

            //- Is the interpolation between consistent meshes? I.e., are all
            //  (non-processor) patches coupled one-to-one?
            bool consistent() const;


        // Source-to-target point finding

            //- Find the target processor and cell associated with a point in a
            //  source cell. See cellsToCells::srcToTgtPoint.
            remote srcToTgtPoint(const label srcCelli, const point& p) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
