/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::zoneGenerators::cellShapes

Description
    A zoneGenerator which constructs a cellZone of the cells
    with the specified shapes

Usage
    \table
        Property     | Description         | Required | Default value
        type         | Type: cellShapes    | yes      |
        name         | Name of the zone    | no       | zoneGenerator name
        shape        | Cell shape          | no       |
        shapes       | List of cell shapes | no       | ()
    \endtable

    To create a cellZone of the hex cells:
    \verbatim
        hexCells
        {
            type        cellShapes;
            shape       hex;
        }
    \endverbatim

    To create a cellZone of the hex and splitHex cells:
    \verbatim
        cells
        {
            type        cellShapes;
            shapes      (hex splitHex);
        }
    \endverbatim

SourceFiles
    cellShapes_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef cellShapes_zoneGenerator_H
#define cellShapes_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellModel;

namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class cellShapes Declaration
\*---------------------------------------------------------------------------*/

class cellShapes
:
    public zoneGenerator
{
    // Private Data

        //- List of cell shape types
        wordList shapeTypes_;

        List<const cellModel*> cellModels_;

        //- Cos of feature angle for polyHedral to be splitHex
        static scalar featureCos;


public:

    //- Runtime type information
    TypeName("cellShapes");


    // Constructors

        //- Construct from dictionary
        cellShapes
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellShapes();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
