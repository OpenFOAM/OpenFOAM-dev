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
    Foam::zoneGenerators::plane

Description
    A zoneGenerator which selects faces based on the adjacent cell centres
    spanning a given plane. The plane is defined by a point and normal vector.

    Additionally, an include entry can be specified. When omitted or set to
    "all", then all faces that meet the criteria are included in the set. When
    set to "closest", just the faces that belong to the closest contiguous
    region to the plane point are included. This latter setting is useful when
    defining face zones through channels on which the flow rate is to be
    computed, as it keeps the set local to a single channel.

    The returned faceZone is oriented and provides the corresponding flipMap.

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: face               | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
        point        | Point on plane           | yes      |
        normal       | Normal to plane          | yes      |
        include      | Include option: all or closest | yes      |
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To create the faceZone \c mid containing the faces based on the adjacent
    cell centres spanning a given plane:
    \verbatim
        mid
        {
            type        plane;
            point       (0.41 0 0);
            normal      (1 0 0);
        }
    \endverbatim

SourceFiles
    planeZoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef planeZoneGenerator_H
#define planeZoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                    Class plane Declaration
\*---------------------------------------------------------------------------*/

class plane
:
    public zoneGenerator
{
public:

    // Public Enumerations

        //- Enumeration for what to include
        enum class include
        {
            all,
            closest
        };

        //- Include option names
        static const NamedEnum<include, 2> includeNames;


private:

    // Private Data

        //- Point on the plane
        const vector point_;

        //- Normal to the plane
        const vector normal_;

        //- Include option
        const include include_;


public:

    //- Runtime type information
    TypeName("plane");


    // Constructors

        //- Construct from dictionary
        plane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~plane();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
