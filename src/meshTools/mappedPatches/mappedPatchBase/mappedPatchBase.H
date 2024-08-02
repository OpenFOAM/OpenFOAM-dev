/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::mappedPatchBase

Description
    Engine and base class for poly patches which provides interpolative mapping
    between two globally conforming poly patches

    Example:
    \verbatim
        // (Settings inherited from mappedPatchBaseBase to identify the
        // neighbouring region and patch)
        ...

        // Patch-to-patch mapping method. How the interpolation addressing and
        // weights are constructed. Likely to be "matching" if the patches are
        // conformal, or "intersection" if they are not. Optional. If omitted,
        // then it will be assumed that the patches are conformal, and an
        // octree-based method will be used to construct the connections.
        method matching; // nearest, intersection, inverseDistance

        // Tolerance used to check the patches' global conformance (i.e., the
        // conformance of the patches' entire surfaces, not of each pair of
        // faces). Optional. Defaults to a small value.
        matchTol 1e-4;
    \endverbatim

See also
    Foam::mappedPatchBaseBase

SourceFiles
    mappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchBase_H
#define mappedPatchBase_H

#include "mappedPatchBaseBase.H"
#include "patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchBase
:
    public mappedPatchBaseBase
{
protected:

    // Protected data

        //- Are we using the tree mapping method, or a patch-to-patch
        //  intersection engine? The former is typically faster on small,
        //  pairs of patches with identical meshes. The latter is better
        //  parallelised and provides options for connecting patches with
        //  differing surface meshes.
        const bool usingTree_;

        //- Distributor (if using tree)
        mutable autoPtr<distributionMap> treeMapPtr_;

        //- Indices of the neighbouring patch faces who's values
        //  have to be supplied to the distribution map (if using tree)
        mutable labelList treeNbrPatchFaceIndices_;

        //- Patch-to-patch intersection engine (if using patch-to-patch)
        mutable autoPtr<patchToPatch> patchToPatchPtr_;

        //- Default Matching tolerance
        static const scalar defaultMatchTol_;

        //- Matching tolerance
        const scalar matchTol_;

        //- Is the mapping engine up to date? Zero if out of date. One if out
        //  of date (or potentially not out of date) because of motion. Two if
        //  up to date.
        mutable label mappingIsValid_;

        //- If the neighbour's mapping engine (if there is one) up to date?
        mutable label nbrMappingIsValid_;


    // Protected Member Functions

        //- Get the face-areas for this patch
        virtual tmp<vectorField> patchFaceAreas() const;

        //- Get the face-centres for this patch
        virtual tmp<pointField> patchFaceCentres() const;

        //- Get the local points for this patch
        virtual tmp<pointField> patchLocalPoints() const;

        //- Get the face-areas for the neighbour patch
        tmp<vectorField> nbrPatchFaceAreas() const;

        //- Get the face-centres for the neighbour patch
        tmp<pointField> nbrPatchFaceCentres() const;

        //- Get the local points for the neighbour patch
        tmp<pointField> nbrPatchLocalPoints() const;

        //- Return whether or not the mapping is valid
        bool mappingIsValid() const;

        //- Calculate mapping
        void calcMapping() const;

        //- Clear out data on mesh change
        virtual void clearOut() const;

        //- Return whether or not the neighbour patch is of mapped type
        inline bool nbrPatchIsMapped() const;

        //- Get the mapped neighbouring patch
        inline const mappedPatchBase& nbrMappedPatch() const;

        //- Determine whether the neighbouring patch maps from this patch in
        //  the same way as this patch maps from its neighbour
        inline bool symmetric() const;


public:

    //- Runtime type information
    TypeName("mappedPatchBase");


    // Constructors

        //- Construct from patch
        mappedPatchBase(const polyPatch&);

        //- Construct from components
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& nbrRegionName,
            const word& nbrPatchName,
            const cyclicTransform& transform
        );

        //- Construct from dictionary
        mappedPatchBase
        (
            const polyPatch& pp,
            const dictionary& dict,
            const bool defaultTransformIsNone
        );

        //- Construct as copy, resetting patch
        mappedPatchBase(const polyPatch&, const mappedPatchBase&);


    //- Destructor
    virtual ~mappedPatchBase();


    // Member Functions

        // Casting

            //- Cast the given polyPatch to a mappedPatchBase. Handle errors.
            static const mappedPatchBase& getMap(const polyPatch& patch);


        // Edit

            //- Clear out data on mesh change
            void clearOut(const bool move);


        // Distribute

            //- Map/interpolate the neighbour patch field to this patch
            template<class Type>
            tmp<Field<Type>> fromNeighbour(const Field<Type>& nbrFld) const;

            //- Map/interpolate the neighbour patch field to this patch
            template<class Type>
            tmp<Field<Type>> fromNeighbour
            (
                const tmp<Field<Type>>& nbrFld
            ) const;

            //- Map/interpolate patch field to the neighbour patch
            template<class Type>
            tmp<Field<Type>> toNeighbour(const Field<Type>& fld) const;

            //- Map/interpolate patch field to the neighbour patch
            template<class Type>
            tmp<Field<Type>> toNeighbour
            (
                const tmp<Field<Type>>& fld
            ) const;


        // I/O

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
