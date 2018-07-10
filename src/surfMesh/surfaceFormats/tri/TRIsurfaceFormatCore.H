/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fileFormats::TRIsurfaceFormatCore

Description
    Internal class used by the TRIsurfaceFormat

SourceFiles
    TRIsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef TRIsurfaceFormatCore_H
#define TRIsurfaceFormatCore_H

#include "surfaceFormatsCore.H"
#include "triFace.H"

#include "IFstream.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class TRIsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class TRIsurfaceFormatCore
:
    public surfaceFormatsCore
{
    // Private Data

        bool sorted_;

        //- The points supporting the facets
        pointField points_;

        //- The zones associated with the faces
        List<label> zoneIds_;

        //- The solid count, in the order of their first appearance
        List<label> sizes_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        TRIsurfaceFormatCore(const TRIsurfaceFormatCore&);

        //- Disallow default bitwise assignment
        void operator=(const TRIsurfaceFormatCore&);

        bool read(const fileName&);


public:

    // Constructors

        //- Read from file, filling in the information
        TRIsurfaceFormatCore(const fileName&);


    //- Destructor
    ~TRIsurfaceFormatCore();


    // Member Functions

        //- File read was already sorted
        bool sorted() const
        {
            return sorted_;
        }

        //- Flush all values
        void clear()
        {
            sorted_ = true;
            points_.clear();
            zoneIds_.clear();
            sizes_.clear();
        }

        //- Return full access to the points
        pointField& points()
        {
            return points_;
        }

        //- Return full access to the zones
        List<label>& zoneIds()
        {
            return zoneIds_;
        }

        //- The list of zone sizes in the order of their first appearance
        List<label>& sizes()
        {
            return sizes_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
