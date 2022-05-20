/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::extendedCellToFaceStencil

Description
    Calculates/contains the extended cell-to-face stencil.

    The stencil is a list of indices into either cells or boundary faces
    in a compact way. (element 0 is owner, 1 is neighbour). The index numbering
    is
    - cells first
    - then all (non-empty patch) boundary faces

    When used in evaluation is a two stage process:
    - collect the data (cell data and non-empty boundaries) into a
    single field
    - (parallel) distribute the field
    - sum the weights*field.

SourceFiles
    extendedCellToFaceStencil.C
    extendedCellToFaceStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCellToFaceStencil_H
#define extendedCellToFaceStencil_H

#include "distributionMap.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class extendedCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCellToFaceStencil
{
protected:

    // Protected data

        const polyMesh& mesh_;


    // Protected member functions

        //- Write some statistics about stencil
        static void writeStencilStats
        (
            Ostream& os,
            const labelListList& stencil,
            const distributionMap& map
        );


public:

    // Declare name of the class and its debug switch
    ClassName("extendedCellToFaceStencil");


    // Constructors

        //- Construct from mesh
        explicit extendedCellToFaceStencil(const polyMesh&);


    // Member Functions

        //- Use map to get the data into stencil order
        template<class T>
        static void collectData
        (
            const distributionMap& map,
            const labelListList& stencil,
            const GeometricField<T, fvPatchField, volMesh>& fld,
            List<List<T>>& stencilFld
        );

        //- Sum vol field contributions to create face values
        template<class Type>
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        weightedSum
        (
            const distributionMap& map,
            const labelListList& stencil,
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<scalar>>& stencilWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extendedCellToFaceStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
