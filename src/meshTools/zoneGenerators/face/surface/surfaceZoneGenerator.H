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
    Foam::zoneGenerators::surfaceZoneGenerator

Description
    A faceZone zoneGenerator which selects faces based on the intersection of
    a surface with the vector between the adjacent cell centres.

    The returned faceZone is unoriented and provided the corresponding flipMap.

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: face               | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
        surface      | searchableSurface type   | yes      |
        tol          | cell-cell vector intersection tolerance | no | rootSmall
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To create the faceZone \c baffles containing the faces based on the adjacent
    cell centres spanning the given triangulated surface:
    \verbatim
        baffles
        {
            type        surface;
            surface     triSurface;
            file        "baffles.obj";
        }
    \endverbatim

    To create the faceZone \c rotatingZone containing the faces based on the
    adjacent cell centres spanning the given cylinder:
    \verbatim
        rotatingZone
        {
            type        surface;
            surface     searchableCylinder;
            point1      (0 0 0.54);
            point2      (0 0 0.66);
            radius      0.2;
        }
    \endverbatim

See also
    searchableSurface
    zoneGenerators::planeZoneGenerator

SourceFiles
    surfaceZoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceZoneGenerator_H
#define surfaceZoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;

namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                    Class surfaceZoneGenerator Declaration
\*---------------------------------------------------------------------------*/

class surfaceZoneGenerator
:
    public zoneGenerator
{
    // Private Data

        //- Surface
        autoPtr<searchableSurface> surfacePtr_;

        //- Tolerance with which to extend the cell-cell vectors
        const scalar tol_;


public:

    //- Runtime type information
    TypeName("surface");


    // Constructors

        //- Construct from dictionary
        surfaceZoneGenerator
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceZoneGenerator();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
