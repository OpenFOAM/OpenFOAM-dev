/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::fileFormats::STLsurfaceFormat

Description
    Provide a means of reading/writing STL files (ASCII and binary).

Note
    For efficiency, the zones are sorted before creating the faces.
    The class is thus derived from MeshedSurface.

SourceFiles
    STLsurfaceFormat.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLsurfaceFormat_H
#define STLsurfaceFormat_H

#include "STLsurfaceFormatCore.H"
#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class STLsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STLsurfaceFormat
:
    public MeshedSurface<Face>
{
    // Private Member Functions

        //- Write Face (ASCII)
        static inline void writeShell
        (
            Ostream&,
            const pointField&,
            const Face&
        );

        //- Write Face (BINARY)
        static inline void writeShell
        (
            ostream&,
            const pointField&,
            const Face&,
            const label zoneI
        );

        //- Disallow default bitwise copy construct
        STLsurfaceFormat(const STLsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const STLsurfaceFormat<Face>&);


public:

    // Constructors

        //- Construct from file name
        STLsurfaceFormat(const fileName&);


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face> > New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face> >
            (
                new STLsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~STLsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy (as ASCII)
        static void writeAscii
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write surface mesh components by proxy (as BINARY)
        static void writeBinary
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write surface mesh components by proxy
        //  as ASCII or BINARY, depending on the extension
        static void write(const fileName&, const MeshedSurfaceProxy<Face>&);

        //- Write UnsortedMeshedSurface (as ASCII) sorted by zone
        static void writeAscii
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface (as BINARY) unsorted by zone
        static void writeBinary
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  as ASCII or BINARY, depending on the extension
        static void write(const fileName&, const UnsortedMeshedSurface<Face>&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "STLsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
