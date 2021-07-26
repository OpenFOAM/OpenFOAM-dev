/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::regionToFace

Description
    A topoSetSource to select faces belonging to topological connected region
    (that contains given point)

SourceFiles
    regionToFace.C

\*---------------------------------------------------------------------------*/

#ifndef regionToFace_H
#define regionToFace_H

#include "topoSetSource.H"
#include "PackedBoolList.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class regionToFace Declaration
\*---------------------------------------------------------------------------*/

class regionToFace
:
    public topoSetSource
{
    // Private Data


        //- Name of set to use
        word setName_;

        //- Coordinate that is nearest/on connected region
        point nearPoint_;

    // Private Member Functions

        //- Walk edge-face-edge
        void markZone
        (
            const indirectPrimitivePatch& patch,
            const label proci,
            const label facei,
            const label zoneI,
            labelList& faceZone
        ) const;

        void combine(topoSet& set, const bool add) const;

public:

    //- Runtime type information
    TypeName("regionToFace");

    // Constructors

        //- Construct from components
        regionToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const point& nearPoint
        );

        //- Construct from dictionary
        regionToFace
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~regionToFace();


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
