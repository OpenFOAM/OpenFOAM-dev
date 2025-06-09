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
    Foam::zoneGenerators::orientZoneGenerator

Description
    A zoneGenerator which sets the face orientation flipMap
    based on the face normal alignment with the given normal.

    The returned faceZone is unoriented and provided the corresponding flipMap.

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: face               | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
        normal       | Reference unit vector    | yes      |
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To orient the faceZone \c cut corresponding to the positive y-direction:
    \verbatim
        orientedCut
        {
            type        orient;
            normal      (0 1 0);
            cut;
        }
    \endverbatim

SourceFiles
    orientZoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef orientZoneGenerator_H
#define orientZoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class orientZoneGenerator Declaration
\*---------------------------------------------------------------------------*/

class orientZoneGenerator
:
    public zoneGenerator
{
    // Private Data

        //- Zone generator supplying the faces
        autoPtr<zoneGenerator> zoneGenerator_;

        //- Reference unit vector face orientations are compared to
        vector normal_;


public:

    //- Runtime type information
    TypeName("orient");


    // Constructors

        //- Construct from dictionary
        orientZoneGenerator
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~orientZoneGenerator();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
