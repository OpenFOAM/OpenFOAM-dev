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
    Foam::zoneGenerators::orient

Description
    A zoneGenerator which sets the face orientation flipMap

    based either on the face normal alignment with the given normal or
    orientation with respect to a given point outside the surface.

    The returned faceZone is oriented and provides the corresponding flipMap.

Usage
    \table
        Property     | Description               | Required | Default value
        type         | Type: face                | yes      |
        name         | Name of the zone          | no       | zoneGenerator name
        normal       | Reference unit vector     | no       |
        outsidePoint | Point outside the surface | no       |
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

    To orient the faceZone \c cut with respect to the point (0 100 0) outside
    the surface:
    \verbatim
        orientedCut
        {
            type        orient;
            point       (0 100 0);
            cut;
        }
    \endverbatim

SourceFiles
    orient_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef orient_zoneGenerator_H
#define orient_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                        Class orient Declaration
\*---------------------------------------------------------------------------*/

class orient
:
    public zoneGenerator
{
    // Private Data

        //- Zone generator supplying the faces
        autoPtr<zoneGenerator> zoneGenerator_;

        //- Reference unit vector face orientations are compared to
        vector normal_;

        //- Seed outside point
        vector outsidePoint_;


    // Private Member Functions

        //- Orientation based on the face normal alignment with the given normal
        boolList normalOrientation(const labelList& faceIndices) const;

        //- Orientation based on the given outside point
        boolList pointOrientation(const faceZone& fZone) const;

        //- Select orientation algorithm based on the user inputs
        boolList orientation(const faceZone& fZone) const;


public:

    //- Runtime type information
    TypeName("orient");


    // Constructors

        //- Construct from dictionary
        orient
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~orient();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
