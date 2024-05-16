/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::searchableSurfaceToFaceZone

Description
    A topoSetSource to select faces based on intersection (of cell-cell
    vector) with a surface.

SourceFiles
    searchableSurfaceToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceToFaceZone_H
#define searchableSurfaceToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;
class faceZoneSet;

/*---------------------------------------------------------------------------*\
                 Class searchableSurfaceToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceToFaceZone
:
    public topoSetSource
{
    // Private Data

        //- Surface
        autoPtr<searchableSurface> surfacePtr_;

        //- Tolerance with which to extend the cell-cell vectors
        const scalar tol_;


    // Private Member Functions

        void combine(faceZoneSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("searchableSurfaceToFaceZone");


    // Constructors

        //- Construct from dictionary
        searchableSurfaceToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceToFaceZone();


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
