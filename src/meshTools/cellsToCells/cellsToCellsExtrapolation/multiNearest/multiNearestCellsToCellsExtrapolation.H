/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::cellsToCellsExtrapolations::multiNearest

Description
    Extrapolation procedure for cells-to-cells interpolations which constructs
    a value for an unmapped cell as a weighted sum of values from multiple
    nearest mapped cells.

    The faces separating mapped and unmapped regions are separated into
    contiguous sets. Each unmapped cell then finds the nearest mapped cell
    adjacent to each set. Values are then interpolated from these multiple
    nearest cells using an inverse-distance weighting.

    This method is useful when mesh-to-mesh mapping opens up a new connection
    between previously disconnected regions of mesh (e.g., when a valve opens
    in an engine simulation). In this case, this method linearly interpolates
    the values along any newly created channels. This prevents discontinuities
    and minimises problematic transients that result from the mapping
    procedure.

SourceFiles
    multiNearestCellsToCellsExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef multiNearestCellsToCellsExtrapolation_H
#define multiNearestCellsToCellsExtrapolation_H

#include "cellsToCellsExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributionMap;

namespace cellsToCellsExtrapolations
{

/*---------------------------------------------------------------------------*\
                            Class multiNearest Declaration
\*---------------------------------------------------------------------------*/

class multiNearest
:
    public cellsToCellsExtrapolation
{
    // Private Data

        //- For each uncoupled cell, the indices of the local cells which are
        //  extrapolated from
        labelListList uncoupledCellLocalCells_;

        //- For each uncoupled cell, the extrapolation weights
        scalarListList uncoupledCellWeights_;

        //- Map from cells to local extrapolation cells
        autoPtr<distributionMap> extrapolationMapPtr_;


    // Private Member Functions

        //- Extrapolate the given field if necessary
        template<class Type>
        void extrapolateType(Field<Type>& fld) const;


public:

    //- Run-time type information
    TypeName("multiNearest");


    // Constructors

        //- Default construct
        multiNearest();


    //- Destructor
    virtual ~multiNearest();


    // Member Functions

        //- Compute the extrapolation addressing if necessary
        void update
        (
            const polyMesh& mesh,
            const PackedBoolList& cellCoupleds
        );

        //- Extrapolate the given field if necessary
        #define declareExtrapolateType(Type, nullArg) \
            virtual void extrapolate(Field<Type>& fld) const;
        FOR_ALL_FIELD_TYPES(declareExtrapolateType);
        #undef declareExtrapolateType
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellsToCellsExtrapolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
