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
    Foam::zoneGenerators::patchZoneGenerator

Description
    A zoneGenerator which constructs a faceZone from a set of patches

    The returned faceZone is unoriented and provided the corresponding flipMap.

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: face              | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        patch        | Name of patch           | no       |
        patches      | List of patch names     | no       | ()
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To convert patch named \c freeSurface into a faceZone of the same name:
    \verbatim
        freeSurface
        {
            type        patch;
            patch       freeSurface;
        }
    \endverbatim

SourceFiles
    patchZoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef patchZoneGenerator_H
#define patchZoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class patchZoneGenerator Declaration
\*---------------------------------------------------------------------------*/

class patchZoneGenerator
:
    public zoneGenerator
{
    // Private Data

        //- Set of patches to convert into faceZone
        labelHashSet patchSet_;


public:

    //- Runtime type information
    TypeName("patch");


    // Constructors

        //- Construct from dictionary
        patchZoneGenerator
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchZoneGenerator();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
