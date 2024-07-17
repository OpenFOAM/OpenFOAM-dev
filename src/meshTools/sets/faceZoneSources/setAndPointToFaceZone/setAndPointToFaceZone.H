/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::setAndPointToFaceZone

Description
    A topoSetSource to select faces based on usage in a faceSet, where an
    behind/inside point is used to orient the faces.

    The given point must be within the mesh, and the provided face set must
    unambiguously separate cells that are connected to the point from those
    that are not. This means that the set must be closed or must span the
    entire width of the domain, from boundary to boundary.

    Whilst more restrictive than setToAndNormalToFaceZone, this method
    allows sets to be reliably converted into zones even when it is not
    possible to specify a normal direction that applies to all the faces;
    for example a set representing a closed surface, or just one with
    substantial variation in its surface normals.

SourceFiles
    setAndPointToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setAndPointToFaceZone_H
#define setAndPointToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class setAndPointToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setAndPointToFaceZone
:
    public topoSetSource
{
    // Private Data

        //- Name of set to use
        word setName_;

        //- Behind/inside point used to orient the faces
        point point_;


public:

    //- Runtime type information
    TypeName("setAndPointToFaceZone");


    // Constructors

        //- Construct from components
        setAndPointToFaceZone
        (
            const polyMesh& mesh,
            const word& setName,
            const vector& normal
        );

        //- Construct from dictionary
        setAndPointToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~setAndPointToFaceZone();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACEZONESOURCE;
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
