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
    Foam::setsToFaceZone

Description
    A topoSetSource to select faces based on usage in a faceSet and cellSet

SourceFiles
    setsToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setsToFaceZone_H
#define setsToFaceZone_H

#include "topoSetSource.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class setsToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setsToFaceZone
:
    public topoSetSource
{
    // Private Data

        //- Name of set to use
        const word faceSetName_;

        //- Name of set to use
        const word cellSetName_;

        //- Whether cellSet is slave cells or master cells
        const Switch flip_;


public:

    //- Runtime type information
    TypeName("setsToFaceZone");


    // Constructors

        //- Construct from components
        setsToFaceZone
        (
            const polyMesh& mesh,
            const word& faceSetName,
            const word& cellSetName,
            const Switch& flip
        );

        //- Construct from dictionary
        setsToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~setsToFaceZone();


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
