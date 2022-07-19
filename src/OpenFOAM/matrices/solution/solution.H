/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::solution

Description
    Selector class for relaxation factors, solver type and solution.

SourceFiles
    solution.C

\*---------------------------------------------------------------------------*/

#ifndef solution_H
#define solution_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class solution Declaration
\*---------------------------------------------------------------------------*/

class solution
:
    public IOdictionary
{
    // Private Data

        //- Dictionary of temporary fields to cache
        dictionary cache_;

        //- Switch for the caching mechanism
        bool caching_;

        //- Dictionary of relaxation factors for all the fields
        dictionary fieldRelaxDict_;

        //- Dictionary of relaxation factors for all the equations
        dictionary eqnRelaxDict_;

        //- Optional default relaxation factor for all the fields
        scalar fieldRelaxDefault_;

        //- Optional default relaxation factor for all the equations
        scalar eqnRelaxDefault_;

        //- Dictionary of solver parameters for all the fields
        dictionary solvers_;


    // Private Member Functions

        //- Read settings from the dictionary
        void read(const dictionary&);


public:

    //- Update from older solver controls syntax
    //  Usually verbose, since we want to know about the changes
    //  Returns the number of settings changed
    static label upgradeSolverDict(dictionary& dict, const bool verbose=true);

    //- Type information
    ClassName("solution");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        solution
        (
            const objectRegistry& obr,
            const fileName& dictName
        );

        //- Disallow default bitwise copy construction
        solution(const solution&) = delete;


    // Member Functions

        // Access

            //- Return true if the given field should be cached
            bool cache(const word& name) const;

            //- Helper for printing cache message
            template<class FieldType>
            static void cachePrintMessage
            (
                const char* message,
                const word& name,
                const FieldType& vf
            );

            //- Return true if the relaxation factor is given for the field
            bool relaxField(const word& name) const;

            //- Return true if the relaxation factor is given for the equation
            bool relaxEquation(const word& name) const;

            //- Return the relaxation factor for the given field
            scalar fieldRelaxationFactor(const word& name) const;

            //- Return the relaxation factor for the given eqation
            scalar equationRelaxationFactor(const word& name) const;

            //- Return the selected sub-dictionary of solvers if the "select"
            //  keyword is given, otherwise return the complete dictionary
            const dictionary& dict() const;

            //- Return the solver controls dictionary
            const dictionary& solversDict() const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solverDict(const word& name) const;


        // Read

            //- Read the solution dictionary
            bool read();


    // Member Operators

        void operator=(const solution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solutionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
