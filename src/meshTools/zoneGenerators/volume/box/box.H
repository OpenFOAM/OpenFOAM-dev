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
    Foam::zoneGenerators::box

Description
    A zoneGenerator which selects points, cells or faces with centres either
    inside or outside a box or list of boxes.

    By default all the points, cells or faces of the mesh are tested for being
    inside or outside the box or boxes but an optional single zone or list of
    zones may be provided and those points, cells or faces are tested instead.
    This provides an efficient method of hierarchical sub-division of space
    where an initial selection of points cells or faces is refined by selecting
    the sub-set inside or outside the given box or boxes rather than having to
    generate another zone to intersect with.

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: box               | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        zoneType     | Type of zone            | yes      |
        select       | Select either the inside or outside       | no  | inside
        moveUpdate   | Switch to update after mesh motion        | no  | false
        box          | Min and max bounding points of the box    | no  |
        boxes        | List of bounding point pairs of the boxes | no  |
        min          | Min bounding point of the box             | no  |
        max          | Max bounding point of the box             | no  |
        zone         | Optional zone to sub-set                  | no  |
        zones        | Optional list of zones to sub-set         | no  |
    \endtable

    A cellZone named \c box1 containing the cells with centres inside a
    single box defined by its minimum and maximum bounding points is generated
    by
    \verbatim
        box1
        {
            type        box;
            zoneType    cell;

            min     (-0.0075 0 -1);
            max     (-0.003 0.0025 1);
        }
    \endverbatim
    or named \c box1 specifically if the zoneGenerator has to be named
    differently:
    \verbatim
        <zoneGenerator name>
        {
            type        box;
            zoneType    cell;
            name        box1;

            min     (-0.0075 0 -1);
            max     (-0.003 0.0025 1);
        }
    \endverbatim
    Alternatively the box can be specified using \c box entry followed by the
    minimum and maximum bounding points:
    \verbatim
        box2
        {
            type        box;
            zoneType    cell;

            box         (-0.0075 0 -1) (-0.003 0.0025 1);
        }
    \endverbatim
    or a list of boxes can be specified using the \c boxes entry followed by a
    list of minimum and maximum bounding point pairs:
    \verbatim
        boxes
        {
            type        box;
            zoneType    cell;

            boxes
            (
                (-0.0075 0 -1) (-0.003 0.0025 1)
                (-0.005 0.001 -1) (-0.001 0.0025 1)
            );
        }
    \endverbatim
    To select all the cells outside a box set the optional \c select control
    to \c outside:
    \verbatim
        box2
        {
            type        box;
            zoneType    cell;
            select      outside;

            box         (-0.0075 0 -1) (-0.003 0.0025 1);
        }
    \endverbatim

See also
    zoneGenerators::volume
    treeBoundBoxList

SourceFiles
    box.C

\*---------------------------------------------------------------------------*/

#ifndef box_zoneGenerator_H
#define box_zoneGenerator_H

#include "volume.H"
#include "treeBoundBoxList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class box Declaration
\*---------------------------------------------------------------------------*/

class box
:
    public volume
{
    // Private Data

        //- Bounding boxes
        treeBoundBoxList boxes_;


    // Private Member Functions

        friend class volume;

        //- Return true if the boxes contain the given point
        inline bool contains(const point& p) const;


public:

    //- Runtime type information
    TypeName("box");


    // Constructors

        //- Construct from dictionary
        box
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~box();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
