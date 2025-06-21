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
    Foam::zoneSet

Description
    Zone container returned by zoneGenerator::generate

    Optionally contains any combination of point, cell or face zones generated
    by the zoneGenerator::generate function which may simply be the zones
    returned by a lookup from the lists of zones stored on the mesh or generated
    by the particular \c zoneGenerator from a geometric search for example.

    Having \c zoneSet optionally hold zones of more than one type is
    particularly useful for the generation of complex coupled or related zone
    systems, for example rotating regions requiring both the \c cellZone of the
    rotating cells and the faces of the outer moving/sliding boundary.

SourceFiles
    zoneSetI.H
    zoneSet.C

See also
    Foam::zoneGenerator

\*---------------------------------------------------------------------------*/

#ifndef zoneSet_H
#define zoneSet_H

#include "NamedEnum.H"
#include "tmpNrc.H"
#include "pointZone.H"
#include "cellZone.H"
#include "faceZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Enumeration defining the zone types
enum class zoneTypes
{
    point,
    cell,
    face
};

//- Named enumeration defining the zone type names
extern const NamedEnum<zoneTypes, 3> zoneTypesNames;

//- Enumeration defining the zone types with an option for all the types
enum class zoneTypesAll
{
    point,
    cell,
    face,
    all
};

//- Named enumeration defining the zone type names
//  with an option for all the types
extern const NamedEnum<zoneTypesAll, 4> zoneTypesAllNames;


/*---------------------------------------------------------------------------*\
                          Class zoneSet Declaration
\*---------------------------------------------------------------------------*/

struct zoneSet
{
    // Public Data

        //- Temporary pointZone
        tmpNrc<pointZone> pZone;

        //- Temporary cellZone
        tmpNrc<cellZone> cZone;

        //- Temporary faceZone
        tmpNrc<faceZone> fZone;


        //- Enumeration defining the zone types
        enum class zoneTypess
        {
            point,
            cell,
            face
        };


    // Constructors

        //- Construct null
        inline zoneSet();

        //- Construct copy transferring content of temporaries
        inline zoneSet(const zoneSet&);

        //- Construct copy moving content
        inline zoneSet(zoneSet&&) = default;

        //- Construct from pointZone pointer
        inline zoneSet(pointZone* pointZonePtr);

        //- Construct from cellZone pointer
        inline zoneSet(cellZone* cellZonePtr);

        //- Construct from faceZone pointer
        inline zoneSet(faceZone* faceZonePtr);

        //- Construct from pointZone, cellZone and faceZone pointers
        inline zoneSet
        (
            pointZone* pointZonePtr,
            cellZone* cellZonePtr,
            faceZone* faceZonePtr
        );

        //- Construct from pointZone
        inline zoneSet(const pointZone& pointZoneRef);

        //- Construct from cellZone
        inline zoneSet(const cellZone& cellZoneRef);

        //- Construct from faceZone
        inline zoneSet(const faceZone& faceZoneRef);

        //- Construct copy transferring content of temporaries if required
        inline zoneSet(const zoneSet&, bool allowTransfer);

        //- Construct and return a clone with a new name
        inline zoneSet clone(const word& name) const;


    // Member Functions

        //- Return true if any of the zone types have been defined
        inline bool valid() const;

        //- Store the temporary zones and return the stored zoneSet
        zoneSet store() const;

        template<class ZoneType>
        inline const ZoneType& zone() const;


    // Member Operators

        //- Assignment transferring the temporary zones
        void operator=(const zoneSet&);

        //- Assignment moving content
        inline zoneSet& operator=(zoneSet&&) = default;
};


template<>
inline const pointZone& zoneSet::zone<pointZone>() const;

template<>
inline const cellZone& zoneSet::zone<cellZone>() const;

template<>
inline const faceZone& zoneSet::zone<faceZone>() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "zoneSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
