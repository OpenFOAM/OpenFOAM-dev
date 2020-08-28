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
    Foam::patchWriter

Description
    Write patch fields

SourceFiles
    patchWriter.C
    patchWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchWriter_H
#define patchWriter_H

#include "pointMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "vtkMesh.H"
#include "indirectPrimitivePatch.H"
#include "PrimitivePatchInterpolation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                         Class patchWriter Declaration
\*---------------------------------------------------------------------------*/

class patchWriter
{
    const vtkMesh& vMesh_;

    const bool binary_;

    const bool nearCellValue_;

    const fileName fName_;

    const labelList patchIDs_;

    std::ofstream os_;

    label nPoints_;

    label nFaces_;

public:

    // Constructors

        //- Construct from components
        patchWriter
        (
            const vtkMesh&,
            const bool binary,
            const bool nearCellValue,
            const fileName&,
            const labelList& patchIDs
        );


    // Member Functions

        std::ofstream& os()
        {
            return os_;
        }

        label nPoints() const
        {
            return nPoints_;
        }

        label nFaces() const
        {
            return nFaces_;
        }

        //- Write cellIDs
        void writePatchIDs();

        //- Write volFields
        template<class Type>
        void write
        (
            const UPtrList<const GeometricField<Type, fvPatchField, volMesh>>&
        );

        //- Write pointFields
        template<class Type>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, pointPatchField, pointMesh>
            >&
        );

        //- Interpolate and write volFields
        template<class Type>
        void write
        (
            const PrimitivePatchInterpolation<primitivePatch>&,
            const UPtrList<const GeometricField<Type, fvPatchField, volMesh>>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
