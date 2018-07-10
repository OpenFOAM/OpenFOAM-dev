/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::faceSelections::searchableSurfaceSelection

Description
    Selects all (internal or coupled) faces intersecting the searchableSurface.

SourceFiles
    searchableSurfaceSelection.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceSelection_H
#define searchableSurfaceSelection_H

#include "faceSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;

namespace faceSelections
{

/*---------------------------------------------------------------------------*\
                     Class searchableSurfaceSelection Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceSelection
:
    public faceSelection
{
    // Private data

        autoPtr<searchableSurface> surfacePtr_;

public:

    //- Runtime type information
    TypeName("searchableSurface");


    // Constructors

        //- Construct from dictionary
        searchableSurfaceSelection
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Clone
        autoPtr<faceSelection> clone() const
        {
            NotImplemented;
            return autoPtr<faceSelection>(nullptr);
        }


    //- Destructor
    virtual ~searchableSurfaceSelection();


    // Member Functions

        virtual void select(const label zoneID, labelList&, boolList&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace faceSelections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
