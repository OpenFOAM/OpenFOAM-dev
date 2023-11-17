/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::solutionControl

Description
    Solution control class

SourceFiles
    solutionControlI.H
    solutionControl.C

\*---------------------------------------------------------------------------*/

#ifndef solutionControl_H
#define solutionControl_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solutionControl Declaration
\*---------------------------------------------------------------------------*/

class solutionControl
:
    public regIOobject
{
protected:

    // Protected data

        //- Reference to the objectRegistry
        const objectRegistry& registry_;

        //- Name of the algorithm
        const word algorithmName_;

        //- Whitespace the same length as the algorithm name
        const string algorithmSpaceStr_;

        //- Flag set true during the final iteration
        mutable bool finalIter_;


    // Protected Member Functions

        //- Return true if the dependencies have been modified
        virtual bool dependenciesModified() const = 0;

        //- Read controls
        virtual bool read() = 0;

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream&) const;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("solutionControl");


    // Constructors

        //- Construct from the time, the name of the algorithm and a registry
        solutionControl
        (
            const objectRegistry& registry,
            const word& algorithmName
        );


    //- Destructor
    virtual ~solutionControl();


    // Member Functions

        //- Return the name of the algorithm
        inline const word& algorithmName() const;

        //- Return some whitespace the same length as the algorithm name
        inline const char* algorithmSpace() const;

        //- Return the dictionary
        virtual const dictionary& dict() const = 0;

        //- Lookup solutionControl from the objectRegistry and return finalIter
        static bool finalIteration(const objectRegistry& registry);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
