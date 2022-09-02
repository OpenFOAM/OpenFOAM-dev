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
    Foam::mappedInternalPatchBase

Description
    Engine which provides mapping from cells to patch faces

    Example:
    \verbatim
        // The name of the region to map from. Optional. Defaults to the same
        // region as the patch.
        neighbourRegion region0;

        // How to offset the patch face centres to the sampling locations.
        // Optional. This will be inferred if there is a single entry of either
        // "distance" or "offset".
        // - normal              : distance along the patch face normals
        // - direction           : specified offset vector
        offsetMode direction;

        // If offsetMode is normal : The normal distance to offset
        distance 1;

        // If offsetMode is direction : The offset vector
        offset (1 0 0);
    \endverbatim

    Note that patch normals point outward, so if offsetMode is \c normal then a
    negative distance will be required in order to map values from inside the
    region.

SourceFiles
    mappedInternalPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedInternalPatchBase_H
#define mappedInternalPatchBase_H

#include "polyBoundaryMesh.H"
#include "distributionMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mappedInternalPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedInternalPatchBase
{
public:

    // Type enumerations

        //- How to project face centres
        enum offsetMode
        {
            NORMAL,     // offset a distance in the face-normal direction
            DIRECTION   // offset with a specified vector
        };

        static const NamedEnum<offsetMode, 2> offsetModeNames_;


protected:

    // Protected data

        //- Patch to map to
        const polyPatch& patch_;

        //- Name of the region to map from
        word nbrRegionName_;

        //- How to offset the patch face centres to the sampling locations
        offsetMode offsetMode_;

        //- Offset distance
        scalar distance_;

        //- Offset vector
        vector offset_;

        //- Distributor
        mutable autoPtr<distributionMap> mapPtr_;

        //- Indices of the cells who's values have to be supplied to the
        //  distribution map
        mutable labelList cellIndices_;


    // Protected Member Functions

        //- Read the offset mode from a dictionary
        offsetMode readOffsetMode(const dictionary& dict) const;

        //- Calculate mapping
        void calcMapping() const;


public:

    //- Runtime type information
    TypeName("mappedInternalPatchBase");


    // Constructors

        //- Construct from patch
        mappedInternalPatchBase(const polyPatch&);

        //- Construct with no offset
        mappedInternalPatchBase
        (
            const polyPatch& pp,
            const word& nbrRegion
        );

        //- Construct from dictionary
        mappedInternalPatchBase(const polyPatch&, const dictionary&);

        //- Construct as copy, resetting patch
        mappedInternalPatchBase
        (
            const polyPatch&,
            const mappedInternalPatchBase&
        );


    //- Destructor
    virtual ~mappedInternalPatchBase();


    // Member Functions

        // Access

            //- Name of the region to map from
            inline const word& nbrRegionName() const;

            //- Is the neighbour region the same as for this patch?
            inline bool sameRegion() const;

            //- Return reference to the parallel distribution map
            inline const distributionMap& map() const;

            //- Return the indices of the cells who's values have to be
            //  supplied to the distribution map
            inline const labelList& cellIndices() const;

            //- Get the region mesh
            const polyMesh& nbrMesh() const;

            //- Get the sample points
            tmp<pointField> samplePoints() const;


        // Edit

            //- Clear out data on mesh change
            void clearOut();


        // Distribute

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            tmp<Field<Type>> distribute(const Field<Type>& fld) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            tmp<Field<Type>> distribute(const tmp<Field<Type>>& fld) const;


        // I/O

            //- Return whether or not the given dictionary contains a
            //  mappedInternalPatchBase specification
            static bool specified(const dictionary& dict);

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedInternalPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedInternalPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
