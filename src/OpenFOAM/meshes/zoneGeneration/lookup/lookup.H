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
    Foam::zoneGenerators::lookup

Description
    A zoneGenerator which looks-up and returns a zoneSet containing point,
    and/or cell and/or faces zones.

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: lookup            | no       | lookup
        name         | Name of the zone(s)     | no       | zoneGenerator name
        zoneType     | Type of zone            | no       | all the zone types
    \endtable

    A full specification for the zoneGenerators::lookup requires the keyword
    which is the name of the zoneGenerator, the type of the zoneGenerator,
    i.e. \c lookup, the name and type of the zone to be looked-up, e.g.
    \verbatim
        <zoneGenerator name>
        {
            type        lookup;
            name        box;
            zoneType    cell;
        }
    \endverbatim
    However usually the name of the zoneGenerator can be chosen to match the
    name of the zone to be looked-up so the name entry defaults to the
    \c zoneGenerator keyword to make the specification simpler, e.g.
    \verbatim
        box
        {
            type        lookup;
            zoneType    cell;
        }
    \endverbatim
    To further simplify the specification of the \c lookup, given that it will
    be commonly used in a zoneGenerator sequence, \c lookup is the default
    type, so the box cellZone can be looked-up with the following simple
    specification:
    \verbatim
        box
        {
            zoneType    cell;
        }
    \endverbatim
    If the \c zoneType is omitted the zoneSet returned by \c lookup will
    contain each of the three zone types corresponding to the zone name that are
    present in the lists of zones.  Usually there will only be the one zone type
    of the given name so the \c zoneType entry can often be omitted and this
    form of default \lookup can simply be specified as the name of
    the zone only, e.g.
    \verbatim
        box;
    \endverbatim
    There is a special \c all zone of each type generated and cached
    automatically by the corresponding zoneList which holds an identity map of
    all the points, cells or faces and is returned by \c lookup using the
    reserved zone name "all", e.g. to return the cellZone corresponding to all
    the cells in the mesh use:
    \verbatim
        all
        {
            zoneType    cell;
        }
    \endverbatim

See also
    zoneGenerator
    zoneSet
    zoneGeneratorList

SourceFiles
    lookup.C

\*---------------------------------------------------------------------------*/

#ifndef lookup_H
#define lookup_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class lookup Declaration
\*---------------------------------------------------------------------------*/

class lookup
:
    public zoneGenerator
{
    // Private Data

        //- The type of zone to lookup (optional)
        zoneTypesAll zoneType_;


    // Private Member Functions

        //- Check if the zone with the given zoneName and type exists
        //  and return a reference otherwise error
        template<class ZoneListType>
        const typename ZoneListType::zoneType& lookupZone
        (
            const ZoneListType& zones,
            const word& zoneName
        ) const;


public:

    //- Runtime type information
    TypeName("lookup");


    // Constructors

        //- Construct from dictionary
        lookup
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~lookup();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
