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
    Foam::fileFormats::FTRsurfaceFormat

Description
    Reading of the (now deprecated and infrequently used)
    Foam Trisurface Format.

Deprecated
    Other formats are better. (deprecated Mar 2009)

SourceFiles
    FTRsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef FTRsurfaceFormat_H
#define FTRsurfaceFormat_H

#include "IFstream.H"
#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class FTRsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FTRsurfaceFormat
:
    public UnsortedMeshedSurface<Face>
{
    // Private classes

        //- read compatibility for ftr patch definitions
        class ftrPatch
        {
            //- Name of patch
            word name_;

            //- Type of patch (ignored since it is usually "empty")
            word type_;

        public:

            const word& name() const
            {
                return name_;
            }

            friend Istream& operator>>(Istream& is, ftrPatch& p)
            {
                is >> p.name_ >> p.type_;
                return is;
            }
        };


    // Private Member Functions

        //- Disallow default bitwise copy construct
        FTRsurfaceFormat(const FTRsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const FTRsurfaceFormat<Face>&);


public:

    // Constructors

        //- Construct from file name
        FTRsurfaceFormat(const fileName&);


    // Selectors

        //- Read file and return surface
        static autoPtr<UnsortedMeshedSurface<Face> > New(const fileName& name)
        {
            return autoPtr<UnsortedMeshedSurface<Face> >
            (
                new FTRsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~FTRsurfaceFormat()
    {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FTRsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
