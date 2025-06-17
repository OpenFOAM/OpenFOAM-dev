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
    Foam::zoneGenerators::write

Description
    A zoneGenerator which writes all existing zones from the mesh

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: write             | yes      |
        zoneType     | Type of zone            | no       | all the zone types
    \endtable

    Examples:
    To write all the zones that are changing:
    \verbatim
        writeAllZones
        {
            type        write;
        }
    \endverbatim
    or to write only the cellZones:
    \verbatim
        writeCellZones
        {
            type        write;
            zoneType    cell;
        }
    \endverbatim
    When called from the zoneGenerators file which is the most likely use
    the simpler syntax to write all the zones that are changing can be used:
    \verbatim
        write;
    \endverbatim
    or to write only the cellZones:
    \verbatim
        write
        {
            zoneType    cell;
        }
    \endverbatim

See also
    zoneGenerators::remove
    zoneGenerators::clear

SourceFiles
    write.C

\*---------------------------------------------------------------------------*/

#ifndef write_H
#define write_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class write Declaration
\*---------------------------------------------------------------------------*/

class write
:
    public zoneGenerator
{
    // Private Member Data

        //- The type of zone to write (optional)
        zoneTypesAll zoneType_;

        //- timeIndex of the last write
        mutable label timeIndex_;


    // Private Member Functions

        template<class ZoneListType>
        void writeZoneType(const ZoneListType& zones) const;


public:

    //- Runtime type information
    TypeName("write");


    // Constructors

        //- Construct from dictionary
        write
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~write();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
