/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::MultiRegionRefs

Description
    Class to wrap a UPtrList of of region-associated objects (meshes, solvers,
    domainDecompositions, ...). Access will return a wrapped reference and will
    set the Info prefix to the region name. The prefix will remain until the
    wrapped reference goes out of scope.

\*---------------------------------------------------------------------------*/

#ifndef MultiRegionRefs_H
#define MultiRegionRefs_H

#include "IOstreams.H"
#include "regionName.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Region>
class MultiRegionRefs;


/*---------------------------------------------------------------------------*\
                          Class RegionRef Declaration
\*---------------------------------------------------------------------------*/

template<class Region>
class RegionRef
{
private:

    // Private Data

        //- The multi-region container
        const MultiRegionRefs<Region>& mrr_;

        //- The region index
        label regioni_;

        //- Reference to the region object
        Region& region_;


public:

    // Constructors

        //- Construct form components
        template<class NonConstRegion>
        RegionRef
        (
            const MultiRegionRefs<NonConstRegion>& mrr,
            const label regioni,
            Region& region
        );

        //- Disallow copy construct
        RegionRef(const RegionRef& rp) = delete;

        //- Move construct
        RegionRef(RegionRef&& rp);


    //- Destructor
    ~RegionRef();


    // Member operators

        //- Cast to reference
        operator Region&() const;

        //- Obtain the reference
        Region& operator()() const;
};


/*---------------------------------------------------------------------------*\
                      Class MultiRegionRefs Declaration
\*---------------------------------------------------------------------------*/

template<class Region>
class MultiRegionRefs
{
    // Private Data

        //- List of region pointers
        UPtrList<Region>& regions_;

        //- Previous prefix to restore on destruction
        const word previousPrefix_;

        //- Table to facilitate access by name
        HashTable<label, word> indices_;


    // Private Member Functions

        //- Are we prefixing?
        bool prefixes() const;

        //- Width of the write prefix
        string::size_type prefixWidth() const;


public:

    // Friendship

        //- Declare friendship with region-reference class
        friend class RegionRef<Region>;

        //- Declare friendship with const-region-reference class
        friend class RegionRef<const Region>;


    // Constructors

        //- Construct from a list of region pointers
        MultiRegionRefs(UPtrList<Region>& regions);


    //- Destructor
    ~MultiRegionRefs();


    // Member Functions

        //- Return the size
        label size() const;


    // Member Operators

        //- Const-access a region
        RegionRef<const Region> operator[](const label regioni) const;

        //- Access a region
        RegionRef<Region> operator[](const label regioni);

        //- Const-access a region
        RegionRef<const Region> operator[](const word& regionName) const;

        //- Access a region
        RegionRef<Region> operator[](const word& regionName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MultiRegionRefs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
