/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::tecplotWriter

Description
    Write binary tecplot files using tecio.

SourceFiles
    tecplotWriter.C
    tecplotWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef tecplotWriter_H
#define tecplotWriter_H

#include "TECIO.h"
#include "Time.H"
#include "indirectPrimitivePatch.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class tecplotWriter Declaration
\*---------------------------------------------------------------------------*/

class tecplotWriter
{
    const Time& runTime_;

public:

    // Constructors

        //- Construct from components
        tecplotWriter(const Time&);


    // Member Functions

        void writeInit
        (
            const word& name,
            const string& varNames,
            const fileName&,
            INTEGER4 tecplotFileType
        ) const;

        //- Write mesh as polyhedral zone
        void writePolyhedralZone
        (
            const word& zoneName,
            const INTEGER4 strandID,
            const fvMesh& mesh,
            const List<INTEGER4>& varLocArray,
            INTEGER4 nFaceNodes
        ) const;

        //- Write surface as polygonal zone
        void writePolygonalZone
        (
            const word& zoneName,
            const INTEGER4 strandID,
            const indirectPrimitivePatch& pp,
            const List<INTEGER4>& varLocArray
        ) const;

        //- Write unordered data (or rather 1D ordered)
        void writeOrderedZone
        (
            const word& zoneName,
            INTEGER4 strandID,
            const label n,
            const List<INTEGER4>& varLocArray
        ) const;

        //- Write mesh
        void writeConnectivity(const fvMesh& mesh) const;

        //- Write surface
        void writeConnectivity(const indirectPrimitivePatch& pp) const;

        void writeEnd() const;

        //- Write generic Field
        template<class Type>
        void writeField(const Field<Type>& fld) const;


        //- Get either fvPatchField or patchInternalField
        template<class Type>
        tmp<Field<Type>> getPatchField
        (
            const bool nearCellValue,
            const GeometricField<Type, fvPatchField, volMesh>& vfld,
            const label patchi
        ) const;

        //- Get mixed field: fvsPatchField for boundary faces and
        //  internalField for internal faces.
        template<class Type>
        tmp<Field<Type>> getFaceField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const labelList& faceLabels
        ) const;

        template<class GeoField>
        static wordList getNames(const PtrList<GeoField>&);

        template<class Type>
        static void getTecplotNames
        (
            const wordList& names,
            const INTEGER4 loc,
            string& varNames,
            DynamicList<INTEGER4>& varLocation
        );

        template<class GeoField>
        static void getTecplotNames
        (
            const PtrList<GeoField>& flds,
            const INTEGER4 loc,
            string& varNames,
            DynamicList<INTEGER4>& varLocation
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tecplotWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
