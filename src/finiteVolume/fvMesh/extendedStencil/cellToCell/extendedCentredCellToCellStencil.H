/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::extendedCentredCellToCellStencil

Description

SourceFiles
    extendedCentredCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToCellStencil_H
#define extendedCentredCellToCellStencil_H

#include "extendedCellToCellStencil.H"
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToCellStencil;

/*---------------------------------------------------------------------------*\
              Class extendedCentredCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredCellToCellStencil
:
    public extendedCellToCellStencil
{
    // Private Data

        //- Swap map for getting neighbouring data
        autoPtr<mapDistribute> mapPtr_;

        //- Per cell the stencil.
        labelListList stencil_;


public:

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit extendedCentredCellToCellStencil(const cellToCellStencil&);

        //- Disallow default bitwise copy construction
        extendedCentredCellToCellStencil
        (
            const extendedCentredCellToCellStencil&
        );


    // Member Functions

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            return mapPtr_();
        }

        //- Return reference to the stencil
        const labelListList& stencil() const
        {
            return stencil_;
        }

        //- After removing elements from the stencil adapt the schedule (map).
        void compact();

        //- Use map to get the data into stencil order
        template<class Type>
        void collectData
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            List<List<Type>>& stencilFld
        ) const
        {
            extendedCellToFaceStencil::collectData
            (
                map(),
                stencil(),
                fld,
                stencilFld
            );
        }

        //- Sum vol field contributions to create cell values
        template<class Type, class WeightType>
        tmp
        <
            GeometricField
            <
                typename outerProduct<WeightType, Type>::type,
                fvPatchField,
                volMesh
            >
        > weightedSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<WeightType>>& stencilWeights
        ) const
        {
            return extendedCellToCellStencil::weightedSum
            (
                map(),
                stencil(),
                fld,
                stencilWeights
            );
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredCellToCellStencil&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
