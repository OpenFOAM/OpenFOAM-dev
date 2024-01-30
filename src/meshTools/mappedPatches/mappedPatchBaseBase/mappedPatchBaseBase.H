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
    Foam::mappedPatchBaseBase

Description
    Base class for engines and poly patches which provide mapping between two
    poly patches

    Example:
    \verbatim
        // The name of the region to map from. Optional. Defaults to the same
        // region as the patch.
        neighbourRegion region0;

        // The name of the patch to map from
        neighbourPatch movingWall;

        // Couple group to specify the region and patch to map from. This is an
        // alternative to specifying neighbourRegion and neighbourPatch
        // directly, as shown above.
        coupleGroup baffleGroup;
    \endverbatim

SourceFiles
    mappedPatchBaseBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchBaseBase_H
#define mappedPatchBaseBase_H

#include "coupleGroupIdentifier.H"
#include "cyclicTransform.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mappedPatchBaseBase Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchBaseBase
{
protected:

    // Protected data

        //- Patch to map to
        const polyPatch& patch_;

        //- Couple group for the region/patch to map from
        const coupleGroupIdentifier coupleGroup_;

        //- Name of the region to map from
        mutable word nbrRegionName_;

        //- Name of the patch to map from
        mutable word nbrPatchName_;

        //- The transformation between the patches
        mutable cyclicTransform transform_;


public:

    //- Runtime type information
    TypeName("mappedPatchBaseBase");


    // Constructors

        //- Construct from patch
        mappedPatchBaseBase(const polyPatch&);

        //- Construct from components
        mappedPatchBaseBase
        (
            const polyPatch& pp,
            const word& nbrRegionName,
            const word& nbrPatchName,
            const cyclicTransform& transform
        );

        //- Construct from dictionary
        mappedPatchBaseBase
        (
            const polyPatch& pp,
            const dictionary& dict,
            const bool transformisNone
        );

        //- Construct as copy, resetting patch
        mappedPatchBaseBase(const polyPatch&, const mappedPatchBaseBase&);


    //- Destructor
    virtual ~mappedPatchBaseBase();


    // Member Functions

        // Access

            //- Name of the region to map from
            inline const word& nbrRegionName() const;

            //- Name of the patch to map from
            inline const word& nbrPatchName() const;

            //- The transformation between the patches
            inline const transformer& transform() const;

            //- Is the neighbour region the same as for this patch?
            inline bool sameRegion() const;

            //- Is the neighbour patch the same as this patch?
            inline bool samePatch() const;

            //- Is the neighbour patch the same as this patch with no transform?
            inline bool sameUntransformedPatch() const;

            //- Is the neighbour available?
            bool haveNbr() const;

            //- Get the mesh for the region to map from
            const polyMesh& nbrMesh() const;

            //- Get the patch to map from
            const polyPatch& nbrPolyPatch() const;


        // Casting

            //- Restrict use of the mapper to certain configurations
            struct from
            {
                static const label any = 0;
                static const label sameRegion = 1;
                static const label differentRegion = 2;
                static const label differentPatch = 4;
            };

            //- Validate that the map exists and is appropriate for the given
            //  set of permitted configurations
            template<class PatchField, class FieldType>
            static void validateMapForField
            (
                const PatchField& field,
                const FieldType& iF,
                const dictionary& context,
                const label froms = from::any
            );

            //- Validate that the map is appropriate for the given
            //  set of permitted configurations
            template<class PatchField, class FieldType>
            void validateForField
            (
                const PatchField& field,
                const FieldType& iF,
                const dictionary& context,
                const label froms = from::any
            ) const;


        // I/O

            //- Return whether or not the given dictionary contains a
            //  mappedPatchBaseBase specification
            static bool specified(const dictionary& dict);

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatchBaseBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchBaseBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
