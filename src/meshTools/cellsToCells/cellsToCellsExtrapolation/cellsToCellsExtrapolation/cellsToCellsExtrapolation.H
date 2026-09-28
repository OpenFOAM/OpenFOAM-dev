/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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
    Foam::cellsToCellsExtrapolation

Description
    Base class for cells-to-cells extrapolations. Extrapolations fill in values
    for target cells that were not successfully interpolated because they fall
    outside the bounds of the source mesh.

SourceFiles
    cellsToCellsExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef cellsToCellsExtrapolation_H
#define cellsToCellsExtrapolation_H

#include "polyMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class cellsToCellsExtrapolation Declaration
\*---------------------------------------------------------------------------*/

class cellsToCellsExtrapolation
{
protected:

    // Protected Data

        //- List of uncoupled cells
        labelList uncoupledCells_;

        //- Is extrapolation necessary?
        bool extrapolation_;


    // Protected Member Functions

        //- Update the extrapolation flag
        virtual void update(const PackedBoolList& cellCoupleds);

        //- Get the set of faces that sit between coupled and uncoupled cells
        void getCutFaces
        (
            const polyMesh& mesh,
            const PackedBoolList& cellCoupleds,
            const globalIndex& globalCellIndex,
            labelList& cutFaces,
            labelList& cutFaceCoupledCells
        ) const;


public:

    //- Run-time type information
    TypeName("cellsToCellsExtrapolation");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        cellsToCellsExtrapolation,
        word,
        (),
        ()
    );


    // Constructors

        //- Default construct
        cellsToCellsExtrapolation();


    //- Destructor
    virtual ~cellsToCellsExtrapolation();


    // Selector

        //- Select from name
        static autoPtr<cellsToCellsExtrapolation> New(const word& type);


    // Member Functions

        //- Compute the extrapolation addressing if necessary
        virtual void update
        (
            const polyMesh& mesh,
            const PackedBoolList& cellCoupleds
        ) = 0;

        //- Extrapolate the given field if necessary
        #define declareExtrapolateType(Type, nullArg) \
            virtual void extrapolate(Field<Type>& fld) const = 0;
        FOR_ALL_FIELD_TYPES(declareExtrapolateType);
        #undef declareExtrapolateType
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
