/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::singleCellFvMesh

Description
    fvMesh as subset of other mesh. Consists of one cell and all original
    bounday faces. Useful when manipulating boundary data. Single internal
    cell only needed to be able to manipulate in a standard way.

SourceFiles
    singleCellFvMesh.C
    singleCellFvMeshInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef singleCellFvMesh_H
#define singleCellFvMesh_H

#include "fvPatchFieldMapper.H"
#include "fvMesh.H"
#include "labelListIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class singleCellFvMesh Declaration
\*---------------------------------------------------------------------------*/

class singleCellFvMesh
:
    public fvMesh
{
    // Private data

        const labelListIOList patchFaceAgglomeration_;

        //- From patch faces back to agglomeration or fine mesh
        labelListIOList patchFaceMap_;

        //- From fine mesh faces to coarse mesh
        labelIOList reverseFaceMap_;

        //- From coarse points back to original mesh
        labelIOList pointMap_;

        //- From fine points to coarse mesh
        labelIOList reversePointMap_;


    // Private Member Functions

        //- Calculate agglomerated mesh
        void agglomerateMesh(const fvMesh&, const labelListList&);


        //- Disallow default bitwise copy construct
        singleCellFvMesh(const singleCellFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const singleCellFvMesh&);


public:

        //- Patch field mapper class for agglomerated meshes
        class agglomPatchFieldMapper
        :
            public fvPatchFieldMapper
        {
            // Private data

                const labelListList& addressing_;
                const scalarListList& weights_;
                bool hasUnmapped_;

        public:

                //- Construct given addressing
                agglomPatchFieldMapper
                (
                    const labelListList& addressing,
                    const scalarListList& weights
                )
                :
                    addressing_(addressing),
                    weights_(weights),
                    hasUnmapped_(false)
                {
                    forAll(addressing_, i)
                    {
                        if (addressing_[i].empty())
                        {
                            hasUnmapped_ = true;
                            break;
                        }
                    }
                }

                virtual label size() const
                {
                    return addressing_.size();
                }

                virtual bool direct() const
                {
                    return false;
                }

                bool hasUnmapped() const
                {
                    return hasUnmapped_;
                }

                virtual const labelListList& addressing() const
                {
                    return addressing_;
                }

                virtual const scalarListList& weights() const
                {
                    return weights_;
                }
        };



    // Constructors

        //- Construct from fvMesh and no agglomeration
        singleCellFvMesh(const IOobject& io, const fvMesh&);

        //- Construct from fvMesh and agglomeration of boundary faces.
        //  agglomeration is per patch, per patch face index the agglomeration
        //  the face goes into.
        singleCellFvMesh
        (
            const IOobject& io,
            const fvMesh&,
            const labelListList& patchFaceAgglomeration
        );

        //- Read from IOobject
        singleCellFvMesh(const IOobject& io);

    // Member Functions

        bool agglomerate() const
        {
            return patchFaceAgglomeration_.size() > 0;
        }

        //- From patchFace on this back to original mesh or agglomeration
        const labelListList& patchFaceMap() const
        {
            return patchFaceMap_;
        }

        //- From point on this back to original mesh
        const labelList& pointMap() const
        {
            return pointMap_;
        }

        //- From face on original mesh to face on this
        const labelList& reverseFaceMap() const
        {
            return reverseFaceMap_;
        }

        //- From point on original mesh to point on this (or -1 for removed
        //  points)
        const labelList& reversePointMap() const
        {
            return reversePointMap_;
        }

        //- Map volField. Internal field set to average, patch fields straight
        //  copies.
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleCellFvMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
