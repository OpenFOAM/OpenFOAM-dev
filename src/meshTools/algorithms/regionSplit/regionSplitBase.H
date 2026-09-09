/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::regionSplitBase

Description
    Base class containing shared functionality for mesh-regioning algorithms

SourceFiles
    regionSplitBase.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplitBase_H
#define regionSplitBase_H

#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class regionSplit Declaration
\*---------------------------------------------------------------------------*/

class regionSplitBase
:
    public labelList
{
protected:

    // Protected Data

        //- The number of regions
        label nRegions_;


    // Protected Member Functions

        //- Return a global indexing engine but limited to the local process
        static autoPtr<globalIndex> localGlobalIndex(const label localSize);

        //- Compact the given local region split. Return the number of regions.
        static label compactLocalRegionSplit
        (
            labelList& cellRegion
        );

        //- Compact the given global region split. Return the number of regions.
        static label compactGlobalRegionSplit
        (
            const globalIndex& globalNonCompactRegions,
            labelList& cellRegion
        );


public:

    //- Runtime type information
    ClassName("regionSplitBase");


    // Constructors

        //- Construct given the number of elements
        regionSplitBase(const int size);

        //- Disallow default bitwise copy construction
        regionSplitBase(const regionSplitBase&) = delete;


    // Member Functions

        //- Return total number of regions
        label nRegions() const
        {
            return nRegions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
