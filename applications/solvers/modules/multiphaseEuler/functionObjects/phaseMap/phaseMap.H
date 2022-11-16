/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::functionObjects::phaseMap

Description
    This functionObject writes the phase-fraction map field alpha.map with
    incremental value ranges for each phase
    e.g., with values 0 for water, 1 for air, 2 for oil etc.

    Example of function object specification:
    \verbatim
    phaseMap
    {
        type            phaseMap;
        libs            ("libmultiphaseEulerFoamFunctionObjects.so");
        writeControl    writeTime;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: phaseMap     | yes         |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    phaseMap.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_phaseMap_H
#define functionObjects_phaseMap_H

#include "fvMeshFunctionObject.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class phaseMap Declaration
\*---------------------------------------------------------------------------*/

class phaseMap
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Constant access to phases
        const phaseSystem::phaseModelList& phases_;


public:

    //- Runtime type information
    TypeName("phaseMap");


    // Constructors

        //- Construct from Time and dictionary
        phaseMap
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        phaseMap(const phaseMap&) = delete;


    //- Destructor
    virtual ~phaseMap();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the force fields
        virtual bool execute();

        //- Write the force fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseMap&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
