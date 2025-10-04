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
    Foam::zoneGenerators::patchCells

Description
    A zoneGenerator which constructs a cellZone of the cells adjacent to a set
    of patches

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: patchCells        | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        patch        | Name of patch           | no       |
        patches      | List of patch names     | no       | ()
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To convert patch named \c freeSurface into a cellZone of the same name:
    \verbatim
        freeSurface
        {
            type        patchCells;
            patch       freeSurface;
        }
    \endverbatim

SourceFiles
    patchCells_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef patchCells_zoneGenerator_H
#define patchCells_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class patchCells Declaration
\*---------------------------------------------------------------------------*/

class patchCells
:
    public zoneGenerator
{
    // Private Data

        //- Set of patches to convert into cellZone
        labelHashSet patchSet_;


public:

    //- Runtime type information
    TypeName("patchCells");


    // Constructors

        //- Construct from dictionary
        patchCells
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchCells();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
