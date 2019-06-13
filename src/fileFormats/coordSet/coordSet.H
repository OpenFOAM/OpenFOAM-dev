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
    Foam::coordSet

Description
    Holds list of sampling positions

SourceFiles
    coordSet.C

\*---------------------------------------------------------------------------*/

#ifndef coordSet_H
#define coordSet_H

#include "pointField.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coordSet Declaration
\*---------------------------------------------------------------------------*/

class coordSet
:
    public pointField
{

public:

    // Public data types

        //- Enumeration defining the output format for coordinates
        enum coordFormat
        {
            XYZ,
            X,
            Y,
            Z,
            DISTANCE
        };


private:

        //- String representation of coordFormat enums
        static const NamedEnum<coordFormat, 5> coordFormatNames_;


protected:

        //- Name
        const word name_;

        //- Axis write type
        const coordFormat axis_;

        //- Cumulative distance "distance" write specifier.
        scalarList curveDist_;


public:

    // Constructors

        //- Construct from components
        coordSet
        (
            const word& name,
            const word& axis
        );


        //- Construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            const List<point>& points,
            const scalarList& curveDist
        );


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        word axis() const
        {
            return coordFormatNames_[axis_];
        }

        //- Cumulative distance
        const scalarList& curveDist() const
        {
            return curveDist_;
        }

        //- Is axis specification a vector
        bool hasVectorAxis() const;

        //- Get coordinate of point according to axis specification.
        //  If axis="distance" is the curveDist[index]
        scalar scalarCoord(const label index) const;

        //- Get point according to axis="xyz" specification
        vector vectorCoord(const label index) const;

        Ostream& write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
