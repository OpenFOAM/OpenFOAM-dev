/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::functionObjects::fluidMaxDeltaT

Description
    Returns the maximum time-step evaluated from time-dependent
    maximum Courant number and maximum time-step specifications.

    The \c maxCo and \c maxDeltaT are provided as \c Function1 of time,
    supporting constant, tabulated and functional specifications.

Usage
    Example of a time-dependent maximum Courant number:
    \verbatim
    fluidMaxDeltaT
    {
        type      fluidMaxDeltaT;

        maxCo     table
        (
            (0   0.2)
            (0.6 0.2)
            (0.7 0.1)
        );

        maxDeltaT 1;
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::Function1

SourceFiles
    fluidMaxDeltaT.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMaxDeltaT_H
#define fluidMaxDeltaT_H

#include "fvMeshFunctionObject.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class fluidMaxDeltaT Declaration
\*---------------------------------------------------------------------------*/

class fluidMaxDeltaT
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Max Courant number function/table
        autoPtr<Function1<scalar>> maxCoPtr_;

        //- Max time-step function/table
        autoPtr<Function1<scalar>> maxDeltaTPtr_;


public:

    //- Runtime type information
    TypeName("fluidMaxDeltaT");


    // Constructors

        //- Construct from Time and dictionary
        fluidMaxDeltaT
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fluidMaxDeltaT(const fluidMaxDeltaT&) = delete;


    //- Destructor
    virtual ~fluidMaxDeltaT();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the time step value
        virtual bool write();

        //- Return the maximum time-step for stable operation
        virtual scalar maxDeltaT() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluidMaxDeltaT&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
