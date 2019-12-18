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
    Foam::functionObjects::dsmcFields

Description
    Calculate intensive fields:
    - UMean
    - translationalT
    - internalT
    - overallT
    from averaged extensive fields from a DSMC calculation.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    dsmcFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_dsmcFields_H
#define functionObjects_dsmcFields_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class dsmcFields Declaration
\*---------------------------------------------------------------------------*/

class dsmcFields
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("dsmcFields");


    // Constructors

        //- Construct from Time and dictionary
        dsmcFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        dsmcFields(const dsmcFields&) = delete;


    //- Destructor
    virtual ~dsmcFields();


    // Member Functions

        //- Read the dsmcFields data
        virtual bool read(const dictionary&);

        //- Do not evaluate the state at the start of the run
        virtual bool executeAtStart() const
        {
            return false;
        }

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the DSMC fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dsmcFields&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
