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
    Foam::zoneGenerators::volume

Description
    Abstract zoneGenerator which selects points, cells or faces with centres
    either inside a volume.

    By default all the points, cells or faces of the mesh are tested for being
    inside or outside the volume but an optional single zone or list of zones
    may be provided and those points, cells or faces are tested instead.  This
    provides an efficient method of hierarchical sub-division of space where an
    initial selection of points cells or faces is refined by selecting the
    sub-set inside or outside the given volume rather than having to generate
    another zone to intersect with.

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type of volume          | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        zoneType     | Type of zone            | yes      |
        select       | Select either the inside or outside | no  | inside
        moveUpdate   | Switch to update after mesh motion  | no  | false
        zone         | Optional zone to sub-set            | no  |
        zones        | Optional list of zones to sub-set   | no  |
    \endtable

SourceFiles
    volume.C

\*---------------------------------------------------------------------------*/

#ifndef volume_zoneGenerator_H
#define volume_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class volume Declaration
\*---------------------------------------------------------------------------*/

class volume
:
    public zoneGenerator
{
protected:

    // Protected Member Data

        zoneTypes zoneType_;

        //- Enumeration for inside/outside selection
        enum class selection {inside, outside};

        //- Enumeration names for inside/outside selection
        static const NamedEnum<selection, 2> selectionNames;

        //- Initialisation function for selection enumeration names
        static NamedEnum<selection, 2>::namesType selectionNames_()
        {
            return {"inside", "outside"};
        }

        //- Inside/outside selection
        selection select_;

        //- Optional list of zoneGenerators which provide the zones to sub-set
        //  Default selection is all of the points, cells or faces
        zoneGeneratorList zoneGenerators_;


    // Protected Member Functions

        template<class ZoneType, class UnaryOp, class ZoneGenType>
        labelList select
        (
            const ZoneGenType& zoneGen,
            const vectorField& pts,
            const UnaryOp& uop
        ) const;

        template<class ZoneType, class UnaryOp, class ZoneGenType>
        labelList select
        (
            const ZoneGenType& zoneGen,
            const zoneGeneratorList& zoneGenerators,
            const vectorField& pts,
            const UnaryOp& uop
        ) const;

        template<class ZoneType, class ZoneGenType>
        labelList selectOp
        (
            const ZoneGenType& zoneGen,
            const zoneGeneratorList& zoneGenerators,
            const vectorField& pts
        ) const;

        template<class UnaryOp, class ZoneGenType>
        labelList select
        (
            const ZoneGenType& zoneGen,
            const zoneGeneratorList& zoneGenerators,
            const vectorField& pts,
            boolList& flipMap,
            const UnaryOp& uop
        ) const;

        template<class ZoneGenType>
        labelList selectOp
        (
            const ZoneGenType& zoneGen,
            const zoneGeneratorList& zoneGenerators,
            const vectorField& pts,
            boolList& flipMap
        ) const;

        //- Generate the zone of the given type
        template<class ZoneGenType>
        zoneSet generate(const ZoneGenType& zone) const;


public:

    // Constructors

        //- Construct from dictionary
        volume
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volume();


    // Member Functions

        virtual zoneSet generate() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volumeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
