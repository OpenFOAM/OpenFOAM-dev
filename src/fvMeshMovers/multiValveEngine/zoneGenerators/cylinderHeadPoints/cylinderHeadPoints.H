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
    Foam::zoneGenerators::cylinderHeadPoints

Description
    A zoneGenerator which selects the points within the engine cylinder head

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: point              | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To create the pointZone \c cylinderHead of the points within the cylinder
    head:
    \verbatim
        cylinderHead
        {
            type        cylinderHeadPoints;
        }
    \endverbatim

SourceFiles
    cylinderHeadPoints.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderHeadPoints_H
#define cylinderHeadPoints_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                         Class cylinderHeadPoints Declaration
\*---------------------------------------------------------------------------*/

class cylinderHeadPoints
:
    public zoneGenerator
{

public:

    //- Runtime type information
    TypeName("cylinderHeadPoints");


    // Constructors

        //- Construct from dictionary
        cylinderHeadPoints
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cylinderHeadPoints();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
