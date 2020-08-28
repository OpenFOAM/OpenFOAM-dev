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
    Foam::fileFormats::VTKsurfaceFormat

Description
    Provide a means of reading/writing VTK legacy format.
    The output is never sorted by zone.

SourceFiles
    VTKsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKsurfaceFormat_H
#define VTKsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "VTKsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class VTKsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class VTKsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTKsurfaceFormatCore
{
    // Private Member Functions

        //- Write header information about number of polygon points
        static void writeHeaderPolygons(Ostream&, const UList<Face>&);


public:

    // Constructors

        //- Construct from file name
        VTKsurfaceFormat(const fileName&);

        //- Disallow default bitwise copy construction
        VTKsurfaceFormat(const VTKsurfaceFormat<Face>&) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face>> New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face>>
            (
                new VTKsurfaceFormat<Face>(name)
            );
        }

    //- Destructor
    virtual ~VTKsurfaceFormat()
    {}


    // Member Functions

        // Write

        //- Write surface mesh components by proxy
        static void write(const fileName&, const MeshedSurfaceProxy<Face>&);

        //- Write UnsortedMeshedSurface, the output remains unsorted
        static void write(const fileName&, const UnsortedMeshedSurface<Face>&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object file
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VTKsurfaceFormat<Face>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VTKsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
