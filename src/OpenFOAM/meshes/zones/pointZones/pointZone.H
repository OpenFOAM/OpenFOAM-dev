/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::pointZone

Description
    A subset of mesh points.
    The labels of points in the zone can be obtained from the addressing()
    list.

    For quick check whether a point belongs to the zone use the lookup
    mechanism in pointZoneList, where all the zoned points are registered
    with their zone number.

SourceFiles
    pointZone.C

\*---------------------------------------------------------------------------*/

#ifndef pointZone_H
#define pointZone_H

#include "Zone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointZoneList;

/*---------------------------------------------------------------------------*\
                          Class pointZone Declaration
\*---------------------------------------------------------------------------*/

class pointZone
:
    public Zone<pointZone, pointZoneList>
{

public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char* const labelsName;


    //- Runtime type information
    TypeName("pointZone");


    // Constructors

        using Zone::Zone;

        //- Construct and return a clone, resetting the mesh zones
        virtual autoPtr<pointZone> clone(const pointZoneList& mz) const
        {
            return autoPtr<pointZone>
            (
                new pointZone(*this, *this, mz)
            );
        }

        //- Construct and return a clone, resetting the point list
        //  and mesh zones
        virtual autoPtr<pointZone> clone
        (
            const labelUList& addr,
            const pointZoneList& mz
        ) const
        {
            return autoPtr<pointZone>
            (
                new pointZone(*this, addr, mz)
            );
        }


    // Member Functions

        //- Return the mesh points
        const pointField& meshCentres() const;

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether zone is synchronised across coupled boundaries. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Update zone using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Member Operators

        using Zone::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
