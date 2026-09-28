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
    Foam::cellsToCellsExtrapolations::none

Description
    An unimplemented extrapolation procedure for cells-to-cells interpolations.
    This simply checks that the interpolation is complete and that
    extrapolation is not necessary and raises an error if that is not the case.

SourceFiles
    noneCellsToCellsExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef noneCellsToCellsExtrapolation_H
#define noneCellsToCellsExtrapolation_H

#include "cellsToCellsExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellsExtrapolations
{

/*---------------------------------------------------------------------------*\
                            Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public cellsToCellsExtrapolation
{
public:

    //- Run-time type information
    TypeName("none");


    // Constructors

        //- Default construct
        none();


    //- Destructor
    virtual ~none();


    // Member Functions

        //- Compute the extrapolation addressing if necessary
        void update
        (
            const polyMesh& mesh,
            const PackedBoolList& cellCoupleds
        );

        //- Extrapolate the given field if necessary
        #define declareExtrapolateType(Type, nullArg)                          \
            virtual void extrapolate(Field<Type>& fld) const                   \
            {                                                                  \
                NotImplemented;                                                \
            }
        FOR_ALL_FIELD_TYPES(declareExtrapolateType);
        #undef declareExtrapolateType
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellsToCellsExtrapolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
