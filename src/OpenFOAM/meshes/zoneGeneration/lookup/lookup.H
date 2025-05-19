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
        type         | zoneGenerator type      | no       | lookup
        name         | Name of the zone(s)     | no       | zoneGenerator name
        zoneType     | Type of zone            | no       | all the zones found
    \endtable

    A full specification for the \c zoneGenerators::lookup requires the keyword
    which is the name of the \c zoneGenerator, the type of the \c zoneGenerator,
    i.e. \c lookup, the name and type of the zone to be looked-up, e.g.
    \verbatim
        lookupBox
        {
            type        lookup;
            name        box;
            zoneType    cell;
        }
    \endverbatim
    However usually the name of the \c zoneGenerator can be chosen to match the
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
    be commonly used in a \c zoneGenerator sequence, \c lookup is the default
    type, so the box \c cellZone can be looked-up with the following simple
    specification:
    \verbatim
        box
        {
            zoneType    cell;
        }
    \endverbatim
    If the \c zoneType is omitted the \c zoneSet returned by \c lookup will
    contain each of the three zone types corresponding to the zone name that are
    present in the lists of zones.  Usually there will only be the one zone type
    of the given name so the \c zoneType entry can often be omitted and this
    form of default \lookup can simply be specified as the name of
    the zone only, e.g.
    \verbatim
        box;
    \endverbatim

See also
    Foam::zoneGenerator
    Foam::zoneSet
    Foam::zoneGeneratorList

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

        //- If the zoneType is not set lookup all the zones of the given name
        bool anyZoneType_;

        //- The type of zone to lookup (optional)
        zoneTypes zoneType_;


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
