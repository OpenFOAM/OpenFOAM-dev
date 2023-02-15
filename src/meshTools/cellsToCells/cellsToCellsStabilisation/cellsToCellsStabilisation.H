/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::cellsToCellsStabilisation

Description
    Stabilisation data and routines for cell-to-cell interpolations

SourceFiles
    cellsToCellsStabilisation.C
    cellsToCellsStabilisationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cellsToCellsStabilisation_H
#define cellsToCellsStabilisation_H

#include "cellsToCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cellsToCellsStabilisation Declaration
\*---------------------------------------------------------------------------*/

class cellsToCellsStabilisation
{
    // Private Data

        //- Is stabilisation occurring?
        bool stabilisation_;

        //- For each cell, the local cell used to stabilise the interpolation
        labelList localStabilisationCells_;

        //- Map from cells to local stabilisation cells
        autoPtr<distributionMap> stabilisationMapPtr_;


public:

    //- Run-time type information
    TypeName("cellsToCellsStabilisation");


    // Constructors

        //- Construct null
        cellsToCellsStabilisation();


    //- Destructor
    virtual ~cellsToCellsStabilisation();


    // Member Functions

        //- Compute the stabilisation addressing if necessary
        void update
        (
            const polyMesh& mesh,
            const PackedBoolList& cellCoupleds
        );

        //- Stabilise the given field if necessary
        template<class Type>
        void stabilise(Field<Type>& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellsToCellsStabilisationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
