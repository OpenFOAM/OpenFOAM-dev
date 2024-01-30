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
    Foam::MultiRegionList

Description
    Class which combines a UPtrList or PtrListr of region-associated objects
    (meshes, solvers, domainDecompositions, ...) with the automatic region
    prefixing provided by MultiRegionRefs.

\*---------------------------------------------------------------------------*/

#ifndef MultiRegionList_H
#define MultiRegionList_H

#include "MultiRegionRefs.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class MultiRegionListBase Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class Container, class Region>
class MultiRegionListBase
:
    private Container<Region>,
    public MultiRegionRefs<Region>
{
public:

    // Constructors

        //- Construct from a list of regions
        MultiRegionListBase(Container<Region>& regions, bool reuse)
        :
            Container<Region>(regions, reuse),
            MultiRegionRefs<Region>(static_cast<Container<Region>&>(*this))
        {}

        //- Construct from a list of regions
        MultiRegionListBase(Container<Region>&& regions)
        :
            Container<Region>(std::move(regions)),
            MultiRegionRefs<Region>(static_cast<Container<Region>&>(*this))
        {}


    // Member Functions

        //- Inherit the size method from the references
        using MultiRegionRefs<Region>::size;


    // Member Operators

        //- Inherit the access operator from the references
        using MultiRegionRefs<Region>::operator[];
};


/*---------------------------------------------------------------------------*\
                      Typedef MultiRegionList Declaration
\*---------------------------------------------------------------------------*/

template<class Region>
using MultiRegionList = MultiRegionListBase<PtrList, Region>;


/*---------------------------------------------------------------------------*\
                      Typedef MultiRegionUList Declaration
\*---------------------------------------------------------------------------*/

template<class Region>
using MultiRegionUList = MultiRegionListBase<UPtrList, Region>;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
