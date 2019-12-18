/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::faceSelections::faceZoneSelection

Description
    Select faces from faceZone

SourceFiles
    faceZoneSelection.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneSelection_H
#define faceZoneSelection_H

#include "faceSelection.H"
#include "pointField.H"
#include "boolList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionSplit;

namespace faceSelections
{

/*---------------------------------------------------------------------------*\
                     Class faceZoneSelection Declaration
\*---------------------------------------------------------------------------*/

class faceZoneSelection
:
    public faceSelection
{
    // Private Data

        //- Name of faceZone
        const word zoneName_;


   // Private Member Functions

public:

    //- Runtime type information
    TypeName("faceZone");


    // Constructors

        //- Construct from dictionary
        faceZoneSelection
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
    virtual ~faceZoneSelection();


    // Member Functions

        //- Apply this selector
        virtual void select(const label zoneID, labelList&, boolList&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace faceSelections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
