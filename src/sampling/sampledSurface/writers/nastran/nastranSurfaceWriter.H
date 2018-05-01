/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::nastranSurfaceWriter

Description
    A surface writer for the Nastran file format - both surface mesh and fields

    formatOptions
    {
        nastran
        {
            // From OpenFOAM field name to Nastran field name
            fields ((pMean PLOAD2));
            // Optional scale
            scale 2.0;
            // Optional format
            format free;    // short, long, free
        }
    };

SourceFiles
    nastranSurfaceWriter.C
    nastranSurfaceWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef nastranSurfaceWriter_H
#define nastranSurfaceWriter_H

#include "surfaceWriter.H"
#include "NamedEnum.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class nastranSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class nastranSurfaceWriter
:
    public surfaceWriter
{
public:

    enum writeFormat
    {
        wfShort,
        wfLong,
        wfFree
    };

    static const NamedEnum<writeFormat, 3> writeFormatNames_;


private:

    // Private data

        //- Write option
        writeFormat writeFormat_;

        //- Map of OpenFOAM field name vs nastran field name
        HashTable<word> fieldMap_;

        //- Scale to apply to values (default = 1.0)
        scalar scale_;


    // Private Member Functions

        //- Initialise the output stream format params
        void formatOS(OFstream& os) const;

        //- Write a co-ordinate
        void writeCoord
        (
            const point& p,
            const label pointi,
            OFstream& os
        ) const;

        //- Write a face element (CTRIA3 or CQUAD4)
        void writeFace
        (
            const word& faceType,
            const labelList& facePts,
            label& nFace,
            OFstream& os
        ) const;

        //- Main driver to write the surface mesh geometry
        void writeGeometry
        (
            const pointField& points,
            const faceList& faces,
            List<DynamicList<face>>& decomposedFaces,
            OFstream& os
        ) const;

        //- Write a face-based value
        template<class Type>
        void writeFaceValue
        (
            const word& nasFieldName,
            const Type& value,
            const label EID,
            OFstream& os
        ) const;

        //- Templated write operation
        template<class Type>
        void writeTemplate
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<Type>& values,
            const bool isNodeValues,
            const bool verbose
        ) const;


public:

    //- Runtime type information
    TypeName("nastran");


    // Constructors

        //- Construct null
        nastranSurfaceWriter();

        //- Construct with some output options
        nastranSurfaceWriter(const dictionary& options);


    //- Destructor
    virtual ~nastranSurfaceWriter();


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
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const bool verbose = false
        ) const;

        //- Write scalarField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<scalar>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const;

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<vector>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const;

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<sphericalTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const;

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<symmTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const;

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<tensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nastranSurfaceWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
