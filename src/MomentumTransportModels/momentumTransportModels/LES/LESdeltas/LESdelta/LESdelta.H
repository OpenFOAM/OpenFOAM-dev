/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::LESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    LESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESdelta_H
#define LESdelta_H

#include "momentumTransportModel.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class LESdelta
{
protected:

    // Protected data

        const momentumTransportModel& momentumTransportModel_;

        volScalarField delta_;


public:

    //- Runtime type information
    TypeName("LESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESdelta,
            dictionary,
            (
                const word& name,
                const momentumTransportModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, momentumTransportModel and dictionary
        LESdelta
        (
            const word& name,
            const momentumTransportModel& turbulence
        );

        //- Disallow default bitwise copy construction
        LESdelta(const LESdelta&) = delete;


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary& dict
        );

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable&
        );


    //- Destructor
    virtual ~LESdelta()
    {}


    // Member Functions

        //- Return momentumTransportModel reference
        const momentumTransportModel& turbulence() const
        {
            return momentumTransportModel_;
        }

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        void operator=(const LESdelta&) = delete;

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
