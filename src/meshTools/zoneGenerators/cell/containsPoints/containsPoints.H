/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::zoneGenerators::containsPoints

Description
    A zoneGenerator which selects the cells containing the given points

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: cell               | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
        points       | List of points           | yes      |
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To create the cellZone \c pointCells of the cells containing the given two
    points:
    \verbatim
        pointCells
        {
            type        containsPoints;
            points      ((0 1 0) (0.1 1 0));
        }
    \endverbatim

SourceFiles
    containsPoints.C

\*---------------------------------------------------------------------------*/

#ifndef containsPoints_H
#define containsPoints_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                         Class containsPoints Declaration
\*---------------------------------------------------------------------------*/

class containsPoints
:
    public zoneGenerator
{
    // Private Data

        //- Points to select cells
        List<point> points_;


public:

    //- Runtime type information
    TypeName("containsPoints");


    // Constructors

        //- Construct from dictionary
        containsPoints
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~containsPoints();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
