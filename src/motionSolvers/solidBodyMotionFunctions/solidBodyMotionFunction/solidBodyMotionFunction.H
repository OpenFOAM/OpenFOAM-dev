/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

Namespace
    Foam::solidBodyMotionFunctions

Description
    Namespace for solid-body motions


Class
    Foam::solidBodyMotionFunction

Description
    Base class for defining solid-body motions

SourceFiles
    solidBodyMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionFunction_H
#define solidBodyMotionFunction_H

#include "Time.H"
#include "dictionary.H"
#include "septernion.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class solidBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionFunction
{
protected:

    // Protected data

        const word name_;
        dictionary SBMFCoeffs_;
        const Time& time_;


public:

    //- Runtime type information
    TypeName("solidBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidBodyMotionFunction,
            dictionary,
            (
                const word& name,
                const dictionary& SBMFCoeffs,
                const Time& runTime
            ),
            (name, SBMFCoeffs, runTime)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        solidBodyMotionFunction
        (
            const word& name,
            const dictionary& dict,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        solidBodyMotionFunction(const solidBodyMotionFunction&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const = 0;


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<solidBodyMotionFunction> New
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& name = "solidBodyMotionFunction"
        );


    //- Destructor
    virtual ~solidBodyMotionFunction();


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict) = 0;

        //- Write in dictionary format
        virtual void writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidBodyMotionFunction&) = delete;
};


void writeEntry(Ostream& os, const solidBodyMotionFunction& f);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
