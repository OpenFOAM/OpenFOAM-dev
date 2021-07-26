/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::planeToFaceZone

Description
    A topoSetSource to select faces based on the adjacent cell centres spanning
    a given plane. The plane is defined by a point and normal vector.

    Additionally, an include entry can be specified. When omitted or set to
    "all", then all faces that meet the criteria are included in the set. When
    set to "closest", just the faces that belong to the closest contiguous
    region to the plane point are included. This latter setting is useful when
    defining face zones through channels on which the flow rate is to be
    computed, as it keeps the set local to a single channel.

SourceFiles
    planeToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef planeToFaceZone_H
#define planeToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceZoneSet;

/*---------------------------------------------------------------------------*\
                       Class planeToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class planeToFaceZone
:
    public topoSetSource
{
    public:

        //- Enumeration for what to include
        enum include
        {
            all,
            closest
        };

        //- Included region names
        static const NamedEnum<include, 2> includeNames_;


private:

    // Private Data


        //- Point on the plane
        const vector point_;

        //- Normal to the plane
        const vector normal_;

        //- Included regions
        const include include_;


    // Private Member Functions

        void combine(faceZoneSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("planeToFaceZone");

    // Constructors

        //- Construct from dictionary
        planeToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~planeToFaceZone();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACESETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
