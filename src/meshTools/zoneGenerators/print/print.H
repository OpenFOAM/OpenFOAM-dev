/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::zoneGenerators::print

Description
    Print a summary of one or more zones to the terminal

    All zone types
        The number of elements (points, faces or cells) in the zone
        The bounding box containing the elements
        The bounding sphere containing the elements

    pointZones
        The average position of all the points

    faceZones
        The total face area magnitude
        The centroid of the faces

    cellZones
        The total volume
        The centroid of the cells

Usage
    \table
        Property     | Description                 | Required | Default value
        type         | Type: print                 | yes      |
        zoneType     | Type of zone                | no       | all
        zone         | Name of the zone to print   | no       |
        zones        | List of zone names to print | no       |
        geometry     | Print the geometry?         | no       | true
    \endtable

    Examples:
    Print a summary of a zone named 'fan'
    \verbatim
        printOneZone
        {
            type        print;

            zone        fan;
        }
    \endverbatim

    Print a summary of a zone named 'fan' and a zone named 'baffle'
    \verbatim
        printTwoZones
        {
            type        print;

            zones
            (
                fan
                baffle
            );
        }
    \endverbatim

    Print a summary of all the zones, by ommitting 'zone' and 'zones'
    \verbatim
        printAllZones
        {
            type        print;
        }
    \endverbatim

    In createZonesDict, shorthand notation can be used to print all zones
    \verbatim
        print;
    \endverbatim

SourceFiles
    print.C

\*---------------------------------------------------------------------------*/

#ifndef print_zoneGenerator_H
#define print_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class print Declaration
\*---------------------------------------------------------------------------*/

class print
:
    public zoneGenerator
{
    // Private Data

        //- Type of zone to print
        zoneTypesAll zoneType_;

        //- List of zone names to print
        autoPtr<wordList> zoneNames_;

        //- Print the geometry
        bool geometry_;


public:

    //- Runtime type information
    TypeName("print");


    // Constructors

        //- Construct from dictionary
        print
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~print();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
