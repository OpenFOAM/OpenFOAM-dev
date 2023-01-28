/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionObjects::fvModel

Description
    functionObject to instantiate and execute an fvModel

    With this \c functionObject it is possible to use the \c clouds \c fvModel
    to track particles without introducing sources into the continuous phase
    transport equations, i.e. one-way coupling.  When executed from the \c
    functions solver module the particles are tracked without evolving the
    continuous phase and without drag or other transfer terms there is no
    coupling, i.e. a pure Lagrangian simulation.

    Example of function object specification:
    \verbatim
    clouds
    {
        type            fvModel;

        executeAtStart  false;

        fvModel
        {
            type    clouds;
            libs    ("liblagrangianParcel.so");
        }
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fvModelFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef fvModelFunctionObject_H
#define fvModelFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "fvModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class fvModel Declaration
\*---------------------------------------------------------------------------*/

class fvModel
:
    public fvMeshFunctionObject
{
    // Private member data

        //- Pointer to the fvModel
        autoPtr<Foam::fvModel> fvModelPtr_;


public:

    //- Runtime type information
    TypeName("fvModel");


    // Constructors

        //- Construct from Time and dictionary
        fvModel
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fvModel(const fvModel&) = delete;


    //- Destructor
    virtual ~fvModel();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute the fvModel
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
