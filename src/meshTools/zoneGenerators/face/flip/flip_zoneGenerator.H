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
    Foam::zoneGenerators::flip

Description
    A zoneGenerator which flips the flipMap of the faceZone returned by the
    given zoneGenerator.

    The returned faceZone is oriented and provides the corresponding flipMap.

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: face               | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
    \endtable

    To create the faceZone floatingObject containing the faces of the
    \c floatingObject patch but flipped to point inwards:
    \verbatim
        floatingObjectBottom
        {
            type        flip;

            floatingObject
            {
                type        patch;
                patch       floatingObject;
            }
        }
    \endverbatim

SourceFiles
    flip_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef flip_zoneGenerator_H
#define flip_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                          Class flip Declaration
\*---------------------------------------------------------------------------*/

class flip
:
    public zoneGenerator
{
    // Private Data

        //- Zone generator supplying the faces
        autoPtr<zoneGenerator> zoneGenerator_;


public:

    //- Runtime type information
    TypeName("flip");


    // Constructors

        //- Construct from dictionary
        flip
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~flip();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
