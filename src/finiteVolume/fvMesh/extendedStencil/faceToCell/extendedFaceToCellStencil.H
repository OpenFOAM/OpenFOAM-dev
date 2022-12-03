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
    Foam::extendedFaceToCellStencil

Description
    Note: transformations on coupled patches not supported. Problem is the
    positions of cells reachable through these patches.

SourceFiles
    extendedFaceToCellStencil.C
    extendedFaceToCellStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFaceToCellStencil_H
#define extendedFaceToCellStencil_H

#include "distributionMap.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class extendedFaceToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedFaceToCellStencil
{
protected:

    // Protected data

        const polyMesh& mesh_;

public:

    // Constructors

        //- Construct from mesh
        explicit extendedFaceToCellStencil(const polyMesh&);


    // Member Functions

        //- Use map to get the data into stencil order
        template<class Type>
        static void collectData
        (
            const distributionMap& map,
            const labelListList& stencil,
            const SurfaceField<Type>& fld,
            List<List<Type>>& stencilFld
        );

        //- Sum surface field contributions to create cell values
        template<class Type>
        static tmp<VolField<Type>> weightedSum
        (
            const distributionMap& map,
            const labelListList& stencil,
            const SurfaceField<Type>& fld,
            const List<List<scalar>>& stencilWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extendedFaceToCellStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
