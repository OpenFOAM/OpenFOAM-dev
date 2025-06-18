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
    Foam::zoneGenerators::remove

Description
    A zoneGenerator which removes pointZones and/or cellZone and/or faceZones
    from the mesh

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: clear             | yes      |
        pointZone    | pointZone name          | no       |
        pointZones   | List of pointZone names | no       |
        cellZone     | cellZone name           | no       |
        cellZones    | List of cellZone names  | no       |
        faceZone     | faceZone name           | no       |
        faceZones    | List of faceZone names  | no       |
    \endtable
    Either a single name or list of names of each zone type can be specified
    and all the named zones are removed from the zone lists held by the mesh.

    Example:
    \verbatim
        removeCellZones
        {
            type        remove;
            cellZones   (cone1 coneAnnulus);
        }
    \endverbatim

See also
    zoneGenerators::clear

SourceFiles
    remove.C

\*---------------------------------------------------------------------------*/

#ifndef remove_zoneGenerator_H
#define remove_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class remove Declaration
\*---------------------------------------------------------------------------*/

class remove
:
    public zoneGenerator
{
    // Private Data

        //- List of pointZone names to remove
        wordList pointZoneNames_;

        //- List of cellZone names to remove
        wordList cellZoneNames_;

        //- List of faceZone names to remove
        wordList faceZoneNames_;


public:

    //- Runtime type information
    TypeName("remove");


    // Constructors

        //- Construct from dictionary
        remove
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~remove();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
