/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::TimeFunction1

Description
    Light wrapper around Function1 to provide a mechanism to update time-based
    entries.

SourceFiles
    TimeFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef TimeFunction1_H
#define TimeFunction1_H

#include "Function1.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class Type>
class TimeFunction1;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TimeFunction1<Type>&
);

/*---------------------------------------------------------------------------*\
                        Class TimeFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TimeFunction1
{
protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Name of the data entry
        const word name_;

        //- The function
        autoPtr<Function1<Type>> function_;


public:

    // Constructor

        //- Construct from time, name and dictionary
        TimeFunction1
        (
            const Time& time,
            const word& name,
            const dictionary& dict
        );

        //- Construct null from time and name
        TimeFunction1
        (
            const Time& time,
            const word& name
        );

        //- Copy constructor
        TimeFunction1(const TimeFunction1<Type>& tde);


    //- Destructor
    virtual ~TimeFunction1();


    // Member Functions

        // Access

            //- Reset entry by re-reading from dictionary
            void reset(const dictionary& dict);


        // Evaluation

            //- Return value as a function of scalar x
            virtual Type value(const scalar x) const;

            //- Integrate between two scalars
            virtual Type integral(const scalar x1, const scalar x2) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const TimeFunction1<Type>& de
        );

        //- Write in dictionary format
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TimeFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
