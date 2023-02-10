/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::WallLocationData

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    WallLocationDataI.H

\*---------------------------------------------------------------------------*/

#ifndef WallLocationData_H
#define WallLocationData_H

#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class transformer;

// Forward declaration of friend functions and operators
template<class WallLocation, class Type>
class WallLocationData;
template<class WallLocation, class Type>
Ostream& operator<<(Ostream&, const WallLocationData<WallLocation, Type>&);
template<class WallLocation, class Type>
Istream& operator>>(Istream&, WallLocationData<WallLocation, Type>&);

/*---------------------------------------------------------------------------*\
                       Class WallLocationData Declaration
\*---------------------------------------------------------------------------*/

template<class WallLocation, class Type>
class WallLocationData
:
    public WallLocation
{
    // Private Data

        //- Data at nearest wall center
        Type data_;


protected:

    // Protected Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point, false
        //  otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const WallLocationData<WallLocation, Type>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Public Typedefs

        //- The type of the stored data
        typedef Type dataType;


    // Constructors

        //- Construct null
        inline WallLocationData();

        //- Construct from data and other geometry
        template<class ... Geometry>
        inline WallLocationData(const Type& data, const Geometry& ... geometry);


    // Member Functions

        // Access

            inline const Type& data() const;

            inline Type& data();

            template<class TrackingData>
            inline const Type& data(TrackingData& td) const;


        // Needed by FaceCellWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const transformer& transform,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<< <WallLocation, Type>
        (
            Ostream&,
            const WallLocationData<WallLocation, Type>&
        );
        friend Istream& operator>> <WallLocation, Type>
        (
            Istream&,
            WallLocationData<WallLocation, Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WallLocationDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
