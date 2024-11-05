/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::fvMeshTools

Description
    A collection of tools for operating on an fvMesh.

SourceFiles
    fvMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTools_H
#define fvMeshTools_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fvMeshTools Declaration
\*---------------------------------------------------------------------------*/

class fvMeshTools
{
    // Private Static Member Functions

        //- Set patchFields according to dictionary
        template<class GeoField>
        static void setPatchFields
        (
            typename GeoField::Mesh& mesh,
            const label patchi,
            const dictionary& patchFieldDict
        );


public:

    // Static Member Functions

        //- Add patch. Inserts patch before all processor patches. Returns the
        //  index of the added patch.
        static label addPatch(fvMesh& mesh, const polyPatch& patch);

        //- Complete adding patches
        static void addedPatches(fvMesh& mesh);

        //- Change patchField on registered fields according to dictionary
        static void setPatchFields
        (
            fvMesh& mesh,
            const label patchi,
            const dictionary& patchFieldDict
        );

        //- Reorder and remove trailing patches. If validBoundary call is
        //  parallel synced and all add the same patch with same settings
        static void reorderPatches
        (
            fvMesh&,
            const labelList& oldToNew,
            const label nPatches,
            const bool validBoundary
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
