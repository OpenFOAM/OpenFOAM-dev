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
    Foam::surfaceWriter

Description
    Base class for surface writers

SourceFiles
    surfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriter_H
#define surfaceWriter_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "pointField.H"
#include "faceList.H"
#include "fileName.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceWriter
{
public:

    //- Runtime type information
    TypeName("surfaceWriter");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            word,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            wordDict,
            (
                const dictionary& optDict
            ),
            (optDict)
        );


    // Selectors

        //- Return a reference to the selected surfaceWriter
        static autoPtr<surfaceWriter> New(const word& writeType);

        //- Return a reference to the selected surfaceWriter
        //  Select with extra write option
        static autoPtr<surfaceWriter> New
        (
            const word& writeType,
            const dictionary& writeOptions
        );


    // Constructors

        //- Construct null
        surfaceWriter();


    //- Destructor
    virtual ~surfaceWriter();


    // Member Functions

        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry()
        {
            return false;
        }

        //- Write single surface geometry to file.
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const bool verbose = false
        ) const = 0;

        //- Write scalarField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const word& fieldName,          // name of field
            const Field<scalar>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const = 0;

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const word& fieldName,          // name of field
            const Field<vector>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const = 0;

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const word& fieldName,          // name of field
            const Field<sphericalTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const = 0;

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const word& fieldName,          // name of field
            const Field<symmTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const = 0;

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const word& fieldName,          // name of field
            const Field<tensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
