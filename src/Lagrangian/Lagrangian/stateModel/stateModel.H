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
    Foam::stateModel

Description
    Base class for models with state

SourceFiles
    stateModel.C

\*---------------------------------------------------------------------------*/

#ifndef stateModel_H
#define stateModel_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class stateModel Declaration
\*---------------------------------------------------------------------------*/

class stateModel
{
    // Private Static Member Functions

        //- Return the state IO object
        static IOobject io
        (
            const word& name,
            const objectRegistry& db,
            const bool processor,
            const bool write
        );


    // Private Classes

        //- State writer class
        class stateModelWriter
        :
            public regIOobject
        {
            // Private Member Data

                //- Reference to the model
                const stateModel& model_;

                //- Is this processor-local state?
                const bool processor_;


        public:

            //- Runtime type information
            TypeName("dictionary");


            // Constructors

                //- Construct from a reference to the model
                stateModelWriter
                (
                    const stateModel& model,
                    const bool processor
                );


            // Member Functions

                //- Write data
                virtual bool writeData(Ostream&) const;
        };


public:

    // Constructors

        //- Construct null
        stateModel();

        //- Disallow default bitwise copy construction
        stateModel(const stateModel&) = delete;


    //- Destructor
    virtual ~stateModel();


    // Member Functions

        //- Return the model name
        virtual const word& name() const = 0;

        //- Return the model's database
        virtual const objectRegistry& db() const = 0;

        //- Write state
        bool writeState(const bool write) const;

        //- Write state
        virtual void writeState(Ostream& os) const;

        //- Write processor state
        virtual void writeProcessorState(Ostream& os) const;

        //- Construct and return the state dictionary for reading
        static dictionary stateDict
        (
            const word& name,
            const objectRegistry& db
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stateModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
