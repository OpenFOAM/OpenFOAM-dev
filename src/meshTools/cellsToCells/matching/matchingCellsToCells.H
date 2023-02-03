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
    Foam::cellsToCellss::matching

Description
    Matching, one-to-one, cells-to-cells interpolation class.

SourceFiles
    matchingCellsToCells.C

\*---------------------------------------------------------------------------*/

#ifndef matchingCellsToCells_H
#define matchingCellsToCells_H

#include "cellsToCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellss
{

/*---------------------------------------------------------------------------*\
                        Class matching Declaration
\*---------------------------------------------------------------------------*/

class matching
:
    public cellsToCells
{
private:

    // Private Member Functions

        //- Return the true if cells intersect
        bool intersect
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            const label srcCelli,
            const label tgtCelli
        ) const;

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

        //- Append to list of src mesh seed indices
        void appendToDirectSeeds
        (
            const polyMesh& srcMesh,
            const polyMesh& tgtMesh,
            boolList& mapFlag,
            labelList& srcTgtSeed,
            DynamicList<label>& srcSeeds,
            label& srcSeedI,
            label& tgtSeedI
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
    TypeName("matching");


    // Constructors

        //- Construct null
        matching();


    //- Destructor
    virtual ~matching();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellsToCellss
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
