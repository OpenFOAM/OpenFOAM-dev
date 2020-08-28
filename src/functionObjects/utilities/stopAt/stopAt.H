/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::functionObjects::stopAt

Description
    Abstract base class for stop conditions.

    Currently the following actions are supported:
    - noWriteNow
    - writeNow
    - nextWrite

SourceFiles
    stopAt.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stopAt_H
#define functionObjects_stopAt_H

#include "functionObject.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class stopAt Declaration
\*---------------------------------------------------------------------------*/

class stopAt
:
    public functionObject
{
public:

        //- Enumeration defining the type of action
        enum class actionType
        {
            noWriteNow,    //!< stop immediately without writing data
            writeNow,      //!< write data and stop immediately
            nextWrite      //!< stop the next time data are written
        };


protected:

    // Protected data

        //- Reference to the Time
        const Time& time_;

        //- Action type names
        static const NamedEnum<actionType, 3> actionTypeNames_;

        //- The type of action
        actionType action_;

        //- Set true when the stop action has executed
        bool stopped_;


    // Protected Member Functions

        //- Return true when the stop condition is achieved
        virtual bool condition() const = 0;


public:

    //- Runtime type information
    TypeName("stopAt");


    // Constructors

        //- Construct from Time and dictionary
        stopAt
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stopAt(const stopAt&) = delete;


    //- Destructor
    virtual ~stopAt();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary&);

        //- Execute, check existence of stopAt file and take action
        virtual bool execute();

        //- Execute, check existence of stopAt file and take action
        virtual bool write();

        //- Execute at the final time-loop, used for cleanup
        virtual bool end();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stopAt&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
