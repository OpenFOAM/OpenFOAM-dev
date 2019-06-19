/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::extendedCentredCellToFaceStencil

Description

SourceFiles
    extendedCentredCellToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToFaceStencil_H
#define extendedCentredCellToFaceStencil_H

#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToFaceStencil;

/*---------------------------------------------------------------------------*\
              Class extendedCentredCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredCellToFaceStencil
:
    public extendedCellToFaceStencil
{
    // Private Data

        //- Swap map for getting neighbouring data
        autoPtr<mapDistribute> mapPtr_;

        //- Per face the stencil.
        labelListList stencil_;


public:

    // Constructors

        //- Construct from uncompacted face stencil
        explicit extendedCentredCellToFaceStencil(const cellToFaceStencil&);

        //- Disallow default bitwise copy construction
        extendedCentredCellToFaceStencil
        (
            const extendedCentredCellToFaceStencil&
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
        template<class T>
        void collectData
        (
            const GeometricField<T, fvPatchField, volMesh>& fld,
            List<List<T>>& stencilFld
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

        //- Sum vol field contributions to create face values
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> weightedSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<scalar>>& stencilWeights
        ) const
        {
            return extendedCellToFaceStencil::weightedSum
            (
                map(),
                stencil(),
                fld,
                stencilWeights
            );
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredCellToFaceStencil&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
