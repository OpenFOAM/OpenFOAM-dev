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
    Foam::zoneGenerators::set

Description
    A zoneGenerator which converts a set into the corresponding zone

    The faceZone generated from a face set is unoriented because face sets do
    not provide the orientation flipMap.

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: set               | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        setName      | Name of the set         | no       | zone name
        zoneType     | Type of set/zone        | yes      |
    \endtable

    Create a cellZone named \c setCells containing the cells in the set:
    \verbatim
        setCells
        {
            type        set;
            zoneType    cell;
            setName     boxCells;
        }
    \endverbatim
    of if the zoneGenerator, cellZone and set all have the name \c boxCells
    this simpler specification can be used:
    \verbatim
        boxCells
        {
            type        set;
            zoneType    cell;
        }
    \endverbatim

SourceFiles
    set.C

\*---------------------------------------------------------------------------*/

#ifndef set_H
#define set_H

#include "zoneGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class set Declaration
\*---------------------------------------------------------------------------*/

class set
:
    public zoneGenerator
{
    // Private Data

        //- Type of the set and corresponding zone
        zoneTypes zoneType_;

        //- Name of the set to be converted
        //  Defaults to the name of the zone
        const word setName_;


public:

    //- Runtime type information
    TypeName("set");


    // Constructors

        //- Construct from dictionary
        set
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~set();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
