/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::cellsToCellss::nearest

Description
    Nearest cells-to-cells interpolation class.

SourceFiles
    nearestCellsToCells.C

\*---------------------------------------------------------------------------*/

#ifndef nearestCellsToCells_H
#define nearestCellsToCells_H

#include "cellsToCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellss
{

/*---------------------------------------------------------------------------*\
                           Class nearest Declaration
\*---------------------------------------------------------------------------*/

class nearest
:
    public cellsToCells
{
private:

    // Private Member Functions

        //- Find indices of overlapping cells in src and tgt meshes - returns
        //  true if found a matching pair
        bool findInitialSeeds
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const label startSeedI,
            label& srcSeedI,
            label& tgtSeedI
        ) const;

        //- Calculate the mesh-to-mesh addressing and weights
        scalar calculateAddressing
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Find the nearest cell on mesh2 for cell1 on mesh1
        void findNearestCell
        (
            const polyMesh& mesh1,
            const polyMesh& mesh2,
            const label cell1,
            label& cell2
        ) const;

        //- Set the next cells for the marching front algorithm
        void setNextNearestCells
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            label& startSeedI,
            label& srcCelli,
            label& tgtCelli,
            boolList& mapFlag,
            const labelList& srcCellIDs
        ) const;

        //- Find a source cell mapped to target cell tgtCelli
        label findMappedSrcCell
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            const label tgtCelli,
            const List<DynamicList<label>>& tgtToSrc
        ) const;


protected:

    // Protected Member Functions

        // Intersection

            //- Calculate the addressing and weights
            virtual scalar calculate
            (
                const polyMesh& srcMesh,
                const polyMesh& tgtMesh
            );

            //- Normalise the weights for a given mesh
            virtual void normalise
            (
                const polyMesh& mesh,
                labelListList& localOtherCells,
                scalarListList& weights
            ) const;


public:

    //- Run-time type information
    TypeName("nearest");


    // Constructors

        //- Construct null
        nearest();


    //- Destructor
    virtual ~nearest();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellsToCellss
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
