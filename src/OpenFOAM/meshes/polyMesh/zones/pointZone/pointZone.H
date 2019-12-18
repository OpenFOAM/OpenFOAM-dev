/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    mechanism in pointZoneMesh, where all the zoned points are registered
    with their zone number.

SourceFiles
    pointZone.C
    pointZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef pointZone_H
#define pointZone_H

#include "zone.H"
#include "pointZoneMeshFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class pointZone;
Ostream& operator<<(Ostream&, const pointZone&);


/*---------------------------------------------------------------------------*\
                          Class pointZone Declaration
\*---------------------------------------------------------------------------*/

class pointZone
:
    public zone
{

protected:

    // Protected data

        //- Reference to zone list
        const pointZoneMesh& zoneMesh_;


public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char * const labelsName;


    //- Runtime type information
    TypeName("pointZone");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointZone,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const pointZoneMesh& zm
            ),
            (name, dict, index, zm)
        );


    // Constructors

        //- Construct from components
        pointZone
        (
            const word& name,
            const labelUList& addr,
            const label index,
            const pointZoneMesh&
        );

        //- Construct from components, transferring contents
        pointZone
        (
            const word& name,
            labelList&& addr,
            const label index,
            const pointZoneMesh&
        );

        //- Construct from dictionary
        pointZone
        (
            const word& name,
            const dictionary&,
            const label index,
            const pointZoneMesh&
        );

        //- Construct given the original zone and resetting the
        //  point list and zone mesh information
        pointZone
        (
            const pointZone&,
            const labelUList& addr,
            const label index,
            const pointZoneMesh&
        );

        //- Construct given the original zone, resetting the
        //  face list and zone mesh information
        pointZone
        (
            const pointZone&,
            labelList&& addr,
            const label index,
            const pointZoneMesh&
        );

        //- Disallow default bitwise copy construction
        pointZone(const pointZone&) = delete;


        //- Construct and return a clone, resetting the zone mesh
        virtual autoPtr<pointZone> clone(const pointZoneMesh& zm) const
        {
            return autoPtr<pointZone>
            (
                new pointZone(*this, *this, index(), zm)
            );
        }

        //- Construct and return a clone, resetting the point list
        //  and zone mesh
        virtual autoPtr<pointZone> clone
        (
            const pointZoneMesh& zm,
            const label index,
            const labelUList& addr
        ) const
        {
            return autoPtr<pointZone>
            (
                new pointZone(*this, addr, index, zm)
            );
        }


    // Selectors

        //- Return a pointer to a new point zone
        //  created on freestore from dictionary
        static autoPtr<pointZone> New
        (
            const word& name,
            const dictionary&,
            const label index,
            const pointZoneMesh&
        );


    //- Destructor
    virtual ~pointZone();


    // Member Functions

        //- Return zoneMesh reference
        const pointZoneMesh& zoneMesh() const;

        //- Helper function to re-direct to zone::localID(...)
        label whichPoint(const label globalPointID) const;

        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const;

        //- Check whether zone is synchronised across coupled boundaries. Return
        //  true if in error.
        virtual bool checkParallelSync(const bool report = false) const;

        //- Correct patch after moving points
        virtual void movePoints(const pointField&)
        {}

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Member Operators

        //- Assignment to zone, clearing demand-driven data
        void operator=(const pointZone&);

        //- Move assignment to zone, clearing demand-driven data
        void operator=(pointZone&&);

        //- Assign addressing, clearing demand-driven data
        void operator=(const labelUList&);

        //- Move addressing, clearing demand-driven data
        void operator=(labelList&&);


    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream&, const pointZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
