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
    the lookup mechanism in meshCellZones, where all the zoned cells are
    registered with their zone number.

SourceFiles
    cellZone.C
    cellZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef cellZone_H
#define cellZone_H

#include "zone.H"
#include "meshCellZonesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cellZone Declaration
\*---------------------------------------------------------------------------*/

class cellZone
:
    public zone
{

protected:

    // Protected data

        //- Reference to zone list
        const meshCellZones& meshZones_;


public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char* const labelsName;


    //- Runtime type information
    TypeName("cellZone");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellZone,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const meshCellZones& mz
            ),
            (name, dict, mz)
        );


    // Constructors

        //- Construct from components
        cellZone
        (
            const word& name,
            const labelUList& addr,
            const meshCellZones&
        );

        //- Construct from components, transferring contents
        cellZone
        (
            const word& name,
            labelList&& addr,
            const meshCellZones&
        );

        //- Construct from dictionary
        cellZone
        (
            const word& name,
            const dictionary&,
            const meshCellZones&
        );

        //- Construct given the original zone and resetting the
        //  cell list and mesh zones information
        cellZone
        (
            const cellZone&,
            const labelUList& addr,
            const meshCellZones&
        );

        //- Construct given the original zone, resetting the
        //  cell list and mesh zones information
        cellZone
        (
            const cellZone&,
            labelList&& addr,
            const meshCellZones&
        );

        //- Disallow default bitwise copy construction
        cellZone(const cellZone&) = delete;


        //- Construct and return a clone, resetting the mesh zones
        virtual autoPtr<cellZone> clone(const meshCellZones& mz) const
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
            const meshCellZones& mz
        ) const
        {
            return autoPtr<cellZone>
            (
                new cellZone(*this, addr, mz)
            );
        }


    // Selectors

        //- Return a pointer to a new cell zone
        //  created on freestore from dictionary
        static autoPtr<cellZone> New
        (
            const word& name,
            const dictionary&,
            const meshCellZones&
        );


    //- Destructor
    virtual ~cellZone();


    // Member Functions

        //- Convenient renaming of zone::localIndex(globalIndex)
        label whichCell(const label globalIndex) const;

        //- Return meshZones reference
        const meshCellZones& meshZones() const;

        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const;

        //- Check whether zone is synchronised across coupled boundaries. Return
        //  true if in error.
        virtual bool checkParallelSync(const bool report = false) const
        {
            return false;
        }

        //- Update zone using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Member Operators

        using zone::operator=;

        //- Assignment to zone, clearing demand-driven data
        void operator=(const cellZone&);

        //- Move assignment to zone, clearing demand-driven data
        void operator=(cellZone&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
