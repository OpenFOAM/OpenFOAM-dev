/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

Description
    Class to encapsulate the topology of a point within a triangle intersection

SourceFiles
    triIntersectLocationI.H
    triIntersectLocationIO.C

\*---------------------------------------------------------------------------*/

#ifndef triIntersectLocation_H
#define triIntersectLocation_H

#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

namespace triIntersect
{

// Forward declaration of friend functions and operators
class location;
bool operator==(const location& a, const location& b);
bool operator!=(const location& a, const location& b);
Ostream& operator<<(Ostream&, const location&);

/*---------------------------------------------------------------------------*\
                          Class location Declaration
\*---------------------------------------------------------------------------*/

class location
:
    private labelPair
{
    private:

        // Private Static Member Functions

            //- Return the null index
            inline static label noi();


        // Private Constructors

            //- Construct from a source index and a target index
            inline location(const label srci, const label tgti);


    public:

        // Static Member Functions

            //- Construct a source point location
            inline static location srcPoint(const label srcPi);

            //- Construct a target point location
            inline static location tgtPoint(const label tgtPi);

            //- Construct a source and target point location
            inline static location srcTgtPoint
            (
                const label srcPi,
                const label tgtPi
            );

            //- Construct an intersection location between a source and a
            //  target edge
            inline static location intersection
            (
                const label srcEi,
                const label tgtEi
            );


        // Constructors

            //- Construct null
            inline location();


        // Member Functions

            //- Return whether the location is null
            inline bool isNull() const;

            //- Return whether the location is a source point
            inline bool isSrcPoint() const;

            //- Return whether the location is a target point
            inline bool isTgtPoint() const;

            //- Return whether the location is a source point and not a target
            //  point
            inline bool isSrcNotTgtPoint() const;

            //- Return whether the location is a target point and not a source
            //  point
            inline bool isTgtNotSrcPoint() const;

            //- Return whether the location is a source point and a target point
            inline bool isSrcAndTgtPoint() const;

            //- Return whether the location is an intersection
            inline bool isIntersection() const;

            //- Return the source point index
            inline label srcPointi() const;

            //- Return the target point index
            inline label tgtPointi() const;

            //- Return the source edge index
            inline label srcEdgei() const;

            //- Return the target edge index
            inline label tgtEdgei() const;


        // Friend Operators

            //- Test if locations are the same
            inline friend bool operator==(const location& a, const location& b);

            //- Test if locations are not the same
            inline friend bool operator!=(const location& a, const location& b);


        // IOstream Operators

            //- Output the location to a stream
            friend Ostream& operator<<(Ostream& os, const location& l);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace triIntersect

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline bool contiguous<triIntersect::location>()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triIntersectLocationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
