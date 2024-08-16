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
    Foam::nonConformalMappedPatchBase

Description
    Base class for poly patches which provides non-conformal mapping between
    two potentially non-globally conforming poly patches

SourceFiles
    nonConformalMappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalMappedPatchBase_H
#define nonConformalMappedPatchBase_H

#include "nonConformalPolyPatch.H"
#include "mappedPatchBaseBase.H"
#include "intersectionPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonConformalMappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class nonConformalMappedPatchBase
:
    public mappedPatchBaseBase
{
protected:

    // Protected data

        //- Non-conformal patch to map to
        const nonConformalPolyPatch& patch_;

        //- Does this side own the patch?
        const bool owner_;

        //- Is the intersection engine up to date? Zero if out of date. One if
        //  out of date (or potentially not out of date) because of motion. Two
        //  if up to date.
        mutable label intersectionIsValid_;

        //- Patch-to-patch intersection engine
        mutable patchToPatches::intersection intersection_;


    // Protected Member Functions

        //- Constructor helper for owner
        bool calcOwner() const;

        //- Neighbour patch
        inline const nonConformalMappedPatchBase& nbrMappedPatch() const;


public:

    //- Runtime type information
    TypeName("nonConformalMappedPatchBase");


    // Constructors

        //- Construct from non-conformal patch
        nonConformalMappedPatchBase(const nonConformalPolyPatch& ncPp);

        //- Construct from components
        nonConformalMappedPatchBase
        (
            const nonConformalPolyPatch& ncPp,
            const word& nbrRegionName,
            const word& nbrPatchName,
            const cyclicTransform& transform
        );

        //- Construct from dictionary
        nonConformalMappedPatchBase
        (
            const nonConformalPolyPatch& ncPp,
            const dictionary& dict,
            const transformType tt
        );

        //- Construct as copy, resetting patch
        nonConformalMappedPatchBase
        (
            const nonConformalPolyPatch& ncPp,
            const nonConformalMappedPatchBase&
        );


    //- Destructor
    virtual ~nonConformalMappedPatchBase();


    // Member Functions

        //- Does this side own the patch?
        bool owner() const;

        //- Access the intersection engine
        const patchToPatches::intersection& intersection() const;

        //- Clear out data on mesh change
        void clearOut(const bool move);

        //- Write as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nonConformalMappedPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
