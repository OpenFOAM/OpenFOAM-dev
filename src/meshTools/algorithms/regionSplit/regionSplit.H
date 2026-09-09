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
    Foam::regionSplit

Description
    This class enumerates disconnected regions of a mesh. It inherits a
    label list which it populates on construction with the region index for
    each cell in the mesh.

    Optionally additional blocked faces and/or explicit connections can be
    provided, and a switch can be used to deactivate communication and keep
    regions on different processors separate.

SourceFiles
    regionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplit_H
#define regionSplit_H

#include "regionSplitBase.H"
#include "boolList.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                         Class regionSplit Declaration
\*---------------------------------------------------------------------------*/

class regionSplit
:
    public regionSplitBase
{
    // Private Data

        //- Reference to the polyMesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Calculate region split in non-compact (global) numbering.
        void calcNonCompactRegionSplit
        (
            const globalIndex& globalFaces,
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections,
            labelList& cellRegion,
            const bool global
        ) const;

        //- Calculate local region split. Return the number of regions.
        label calcLocalRegionSplit
        (
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections,
            labelList& cellRegion
        ) const;

        //- Calculate global region split. Return the number of regions.
        label calcGlobalRegionSplit
        (
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections,
            labelList& cellRegion
        ) const;

        //- Calculate global region split. Return the number of regions.
        label calcRegionSplit
        (
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections,
            labelList& cellRegion,
            const bool global
        ) const;


public:

    //- Runtime type information
    ClassName("regionSplit");


    // Constructors

        //- Construct from mesh
        regionSplit
        (
            const polyMesh& mesh,
            const bool global = Pstream::parRun()
        );

        //- Construct from mesh and whether face is blocked
        //  NOTE: blockedFace has to be consistent across coupled faces!
        regionSplit
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const bool global = Pstream::parRun()
        );

        //- Construct from mesh and whether face is blocked. Additional explicit
        //  connections between normal boundary faces.
        //  NOTE: blockedFace has to be consistent across coupled faces!
        regionSplit
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const List<labelPair>&,
            const bool global = Pstream::parRun()
        );


    // Member Functions

        //- Return reference to the polyMesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
