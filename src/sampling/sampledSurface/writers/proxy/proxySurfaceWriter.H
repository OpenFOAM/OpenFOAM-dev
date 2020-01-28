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
    Foam::proxySurfaceWriter

Description
    A surfaceWriter that writes the geometry via the MeshedSurfaceProxy, but
    which does not support any fields.

SourceFiles
    proxySurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef proxySurfaceWriter_H
#define proxySurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class proxySurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class proxySurfaceWriter
:
    public surfaceWriter
{

    // Private Data

        //- The associated file extension
        word ext_;

public:

    //- Runtime type information
    TypeName("proxy");


    // Constructors

        //- Construct for a given extension
        proxySurfaceWriter(const word& ext);


    //- Destructor
    virtual ~proxySurfaceWriter();


    // Member Functions


        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry()
        {
            return true;
        }


        //- Write single surface geometry to file.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces
        ) const;

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
            const bool isNodeValues
        ) const
        {}

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
            const bool isNodeValues
        ) const
        {}

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
            const bool isNodeValues
        ) const
        {}

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
            const bool isNodeValues
        ) const
        {}

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
            const bool isNodeValues
        ) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
