/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::cellZone

Description
    Named list of cell indices representing a sub-set of the mesh

    Used by mesh-manipulation tools, field initialisation and for sources and
    functionObjects that apply to sub-sets of the cells.

See also
    cellZoneList
    Zone
    zoneGenerator

SourceFiles
    cellZone.C

\*---------------------------------------------------------------------------*/

#ifndef cellZone_H
#define cellZone_H

#include "Zone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellZoneList;

/*---------------------------------------------------------------------------*\
                          Class cellZone Declaration
\*---------------------------------------------------------------------------*/

class cellZone
:
    public Zone<cellZone, cellZoneList>
{

public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char* const labelsName;

        //- Runtime type information
        ClassName("cellZone");


    // Constructors

        using Zone::Zone;

        //- Construct and return a clone
        autoPtr<cellZone> clone() const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, name(), *this, zones())
            );
        }

        //- Construct and return a clone with a new name
        autoPtr<cellZone> clone(const word& name) const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, name, *this, zones())
            );
        }

        //- Construct and return a clone, resetting the mesh zones
        autoPtr<cellZone> clone(const cellZoneList& mz) const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, name(), *this, mz)
            );
        }

        //- Construct and return a clone, resetting the cell list
        //  and mesh zones
        autoPtr<cellZone> clone
        (
            const labelUList& addr,
            const cellZoneList& mz
        ) const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, name(), addr, mz)
            );
        }


    // Member Functions

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether zone is synchronised across coupled boundaries. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const
        {
            return false;
        }

        //- Update zone using the given map
        void topoChange(const polyTopoChangeMap& map);

        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        using Zone::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
