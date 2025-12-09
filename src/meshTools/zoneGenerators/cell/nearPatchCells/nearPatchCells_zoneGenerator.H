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
    Foam::zoneGenerators::nearPatchCells

Description
    A zoneGenerator which constructs a cellZone of the cells within a distance
    from a set of patches

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: nearPatchCells    | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        patch        | Name of patch           | no       |
        patches      | List of patch names     | no       | ()
        distance     | Distance from the patch | yes      |
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To select all cells within one centimetre of a patch named \c walls:
    \verbatim
        nearToWalls
        {
            type        nearPatchCells;
            patch       walls;
            distance    1 [cm];
        }
    \endverbatim

SourceFiles
    nearPatchCells_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef nearPatchCells_zoneGenerator_H
#define nearPatchCells_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                         Class nearPatchCells Declaration
\*---------------------------------------------------------------------------*/

class nearPatchCells
:
    public zoneGenerator
{
    // Private Data

        //- Set of patches
        labelHashSet patchSet_;

        //- Distance from the patch to include in the zone
        scalar distance_;


public:

    //- Runtime type information
    TypeName("nearPatchCells");


    // Constructors

        //- Construct from dictionary
        nearPatchCells
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~nearPatchCells();


    // Member Functions

        //- Generate and return the zoneSet
        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
