/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::fvMeshAdder

Description
    Adds two fvMeshes without using any polyMesh morphing.
    Uses polyMeshAdder.

SourceFiles
    fvMeshAdder.C
    fvMeshAdderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshAdder_H
#define fvMeshAdder_H

#include "polyMeshAdder.H"
#include "volFieldsFwd.H"
#include "fvPatchFieldsFwd.H"
#include "fvsPatchFieldsFwd.H"
#include "fieldMapper.H"
#include "DimensionedField.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class IOobject;
class faceCoupleInfo;
class IOobjectList;
class fvMesh;
class volMesh;
class surfaceMesh;
class mapAddedPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class fvMeshAdder Declaration
\*---------------------------------------------------------------------------*/

class fvMeshAdder
:
    public polyMeshAdder
{
    // Private Member Functions

        //- Calculate map from new patch faces to old patch faces. -1 where
        //  could not map.
        static labelList calcPatchMap
        (
            const label oldStart,
            const label oldSize,
            const labelList& oldToNew,
            const polyPatch& newPatch,
            const label unmappedIndex
        );

        //- Update single volField.
        template<class Type>
        static void MapVolField
        (
            const mapAddedPolyMesh& meshMap,

            VolField<Type>& fld,
            const VolField<Type>& fldToAdd
        );

        //- Update single surfaceField.
        template<class Type>
        static void MapSurfaceField
        (
            const mapAddedPolyMesh& meshMap,

            SurfaceField<Type>& fld,
            const SurfaceField<Type>& fldToAdd
        );

        //- Update single pointField.
        template<class Type>
        static void MapPointField
        (
            const pointMesh& mesh,
            const mapAddedPolyMesh& meshMap,
            const labelListList& oldMeshPoints,

            PointField<Type>& fld,
            const PointField<Type>& fldToAdd
        );

        //- Update single dimensionedField.
        template<class Type>
        static void MapDimField
        (
            const mapAddedPolyMesh& meshMap,

            DimensionedField<Type, volMesh>& fld,
            const DimensionedField<Type, volMesh>& fldToAdd
        );

public:

    // Declare name of the class and its debug switch
    ClassName("fvMeshAdder");


    // Member Functions

        //- Inplace add mesh to fvMesh. Maps all stored fields. Returns map.
        static autoPtr<mapAddedPolyMesh> add
        (
            fvMesh& mesh0,
            const fvMesh& mesh1,
            const faceCoupleInfo& coupleInfo,
            const bool validBoundary = true
        );

        //- Map all volFields of Type
        template<class Type>
        static void MapVolFields
        (
            const mapAddedPolyMesh&,
            const fvMesh& mesh,
            const fvMesh& meshToAdd
        );

        //- Map all surfaceFields of Type
        template<class Type>
        static void MapSurfaceFields
        (
            const mapAddedPolyMesh&,
            const fvMesh& mesh,
            const fvMesh& meshToAdd
        );

        //- Map all surfaceFields of Type
        template<class Type>
        static void MapPointFields
        (
            const mapAddedPolyMesh&,
            const pointMesh& mesh,
            const labelListList& oldMeshPoints,
            const objectRegistry& meshToAdd
        );

        //- Map all DimensionedFields of Type
        template<class Type>
        static void MapDimFields
        (
            const mapAddedPolyMesh&,
            const fvMesh& mesh,
            const fvMesh& meshToAdd
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshAdderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
