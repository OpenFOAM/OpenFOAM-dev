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
    Foam::zoneGenerators::normalZoneGenerator

Description
    A zoneGenerator which selects the faces returned by the given
    zoneGenerator that are aligned with the given normal.

    The returned faceZone is oriented and provides the corresponding flipMap.

Usage
    \table
        Property     | Description              | Required | Default value
        type         | Type: face               | yes      |
        name         | Name of the zone         | no       | zoneGenerator name
        normal       | Reference unit vector    | yes      |
        tol          | cos angle tolerance      | yes      |
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    To create the faceZone floatingObject containing the faces of the
    \c floatingObject patch which point in the y-direction:
    \verbatim
        floatingObjectBottom
        {
            type        normal;

            normal      (0 1 0);
            tol         0.01;

            floatingObject
            {
                type        patch;
                patch       floatingObject;
            }
        }
    \endverbatim

SourceFiles
    normalZoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef normalZoneGenerator_H
#define normalZoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class normalZoneGenerator Declaration
\*---------------------------------------------------------------------------*/

class normalZoneGenerator
:
    public zoneGenerator
{
    // Private Data

        //- Zone generator supplying the faces
        autoPtr<zoneGenerator> zoneGenerator_;

        //- Reference unit vector face orientations are compared to
        vector normal_;

        //- Direction comparison tolerance
        //  cos of angle between normal_ and face normal
        const scalar tol_;


public:

    //- Runtime type information
    TypeName("normal");


    // Constructors

        //- Construct from dictionary
        normalZoneGenerator
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~normalZoneGenerator();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
