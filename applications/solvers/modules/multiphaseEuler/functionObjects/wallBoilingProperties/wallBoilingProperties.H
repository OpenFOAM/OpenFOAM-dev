/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::functionObjects::wallBoilingProperties

Description
    This function looks up wall boiling wall functions and collects and writes
    out out the following data:

      - Bubble departure diameter
      - Bubble departure frequency
      - Nucleation site density
      - Effective liquid fraction at the wall
      - Quenching heat flux
      - Evaporative heat flux

    Example of function object specification:
    \verbatim
    writeWallBoilingProperties
    {
        type            wallBoilingProperties;
        functionObjectLibs ( "libmultiphaseEulerFoamFunctionObjects.so" );
        writeControl    writeTime;
        phase           liquid;
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallBoilingProperties     | yes        |
        phase    | phase name                 | yes         | none
    \endtable

SourceFiles
    wallBoilingProperties.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoilingProperties_H
#define wallBoilingProperties_H

#include "fvMeshFunctionObject.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class wallBoilingProperties Declaration
\*---------------------------------------------------------------------------*/

class wallBoilingProperties
:
    public fvMeshFunctionObject
{
    // Private data

        //- Phase model
        const phaseModel& phase_;

        //- Constant access to phaseSystem
        const phaseSystem& fluid_;


public:

    //- Runtime type information
    TypeName("wallBoilingProperties");


    // Constructors

        //- Construct from Time and dictionary
        wallBoilingProperties
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wallBoilingProperties(const wallBoilingProperties&) = delete;


    //- Destructor
    virtual ~wallBoilingProperties();


    // Member Functions

        //- Read the wallBoilingProperties data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the wallBoilingProperties field
        virtual bool execute();

        //- Write the wallBoilingProperties field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wallBoilingProperties&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
