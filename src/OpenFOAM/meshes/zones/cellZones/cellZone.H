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
    Foam::cellZone

Description
    A subset of mesh cells.

    Currently set up as an indirect list but will be extended to use a
    primitive mesh.  For quick check whether a cell belongs to the zone use
    the lookup mechanism in cellZoneList, where all the zoned cells are
    registered with their zone number.

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

    typedef cellZoneList ZonesType;

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char* const labelsName;

    //- Runtime type information
    TypeName("cellZone");


    // Constructors

        using Zone::Zone;

        //- Construct and return a clone, resetting the mesh zones
        virtual autoPtr<cellZone> clone(const cellZoneList& mz) const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, *this, mz)
            );
        }

        //- Construct and return a clone, resetting the cell list
        //  and mesh zones
        virtual autoPtr<cellZone> clone
        (
            const labelUList& addr,
            const cellZoneList& mz
        ) const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, addr, mz)
            );
        }


    // Member Functions

        //- Return the mesh cell centres
        const pointField& meshCentres() const;

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether zone is synchronised across coupled boundaries. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const
        {
            return false;
        }

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
