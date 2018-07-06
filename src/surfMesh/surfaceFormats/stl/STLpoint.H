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
    Foam::STLpoint

Description
    A vertex point representation for STL files.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef STLpoint_H
#define STLpoint_H

#include "point.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class STLpoint Declaration
\*---------------------------------------------------------------------------*/

class STLpoint
:
    public Vector<float>
{

public:

    // Constructors

        //- Construct null
        inline STLpoint()
        {}

        //- Construct from components
        inline STLpoint(float x, float y, float z)
        :
            Vector<float>(x, y, z)
        {}

        //- Construct from point
        inline STLpoint(const point& pt)
        :
            Vector<float>(float(pt.x()), float(pt.y()), float(pt.z()))
        {}

        //- Construct from istream
        inline STLpoint(Istream& is)
        :
            Vector<float>(is)
        {}


    // Member Operators

        #if defined(WM_DP) || defined(WM_LP)
        //- Conversion to double-precision point
        inline operator point() const
        {
            return point(x(), y(), z());
        }
        #endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
