/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::functionObjects::abort

Description
    Watches for presence of the named file in the $FOAM_CASE directory
    and aborts the calculation if it is present.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite

SourceFiles
    abort.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_abort_H
#define functionObjects_abort_H

#include "functionObject.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class abort Declaration
\*---------------------------------------------------------------------------*/

class abort
:
    public functionObject
{
public:

    // Public data

        //- Enumeration defining the type of action
        enum class actionType
        {
            noWriteNow,    //!< stop immediately without writing data
            writeNow,      //!< write data and stop immediately
            nextWrite      //!< stop the next time data are written
        };

private:

    // Private Data

        //- Reference to the Time
        const Time& time_;

        //- The fully-qualified name of the abort file
        fileName abortFile_;

        //- Action type names
        static const NamedEnum<actionType, 3> actionTypeNames_;

        //- The type of action
        actionType action_;


    // Private Member Functions

        //- Remove abort file.
        void removeFile() const;


public:

    //- Runtime type information
    TypeName("abort");


    // Constructors

        //- Construct from Time and dictionary
        abort
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        abort(const abort&) = delete;


    //- Destructor
    virtual ~abort();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary&);

        //- Execute, check existence of abort file and take action
        virtual bool execute();

        //- Execute, check existence of abort file and take action
        virtual bool write();

        //- Execute at the final time-loop, used for cleanup
        virtual bool end();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const abort&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
