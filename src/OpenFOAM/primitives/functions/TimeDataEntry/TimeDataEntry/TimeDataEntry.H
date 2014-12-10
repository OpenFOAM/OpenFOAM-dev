/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::TimeDataEntry

Description
    Light wrapper around DataEntry to provide a mechanism to update time-based
    entries.

SourceFiles
    TimeDataEntry.C

\*---------------------------------------------------------------------------*/

#ifndef TimeDataEntry_H
#define TimeDataEntry_H

#include "DataEntry.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class Type>
class TimeDataEntry;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TimeDataEntry<Type>&
);

/*---------------------------------------------------------------------------*\
                       Class TimeDataEntry Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TimeDataEntry
{

protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Name of the data entry
        const word name_;

        //- The underlying DataEntry
        autoPtr<DataEntry<Type> > entry_;


public:

    // Constructor

        //- Construct from entry name
        TimeDataEntry
        (
            const Time& t,
            const word& name,
            const dictionary& dict
        );

        //- Construct null from entry name
        TimeDataEntry
        (
            const Time& t,
            const word& entryName
        );

        //- Copy constructor
        TimeDataEntry(const TimeDataEntry<Type>& tde);


    //- Destructor
    virtual ~TimeDataEntry();


    // Member Functions

        // Access

            //- Reset entry by re-reading from dictionary
            void reset(const dictionary& dict);

            //- Return the name of the entry
            const word& name() const;


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const scalar x) const;

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const TimeDataEntry<Type>& de
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TimeDataEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
