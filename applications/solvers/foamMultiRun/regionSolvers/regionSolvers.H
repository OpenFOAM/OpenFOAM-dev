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
    Foam::regionSolvers

Description
    Class to hold the lists of region meshes and solvers

    Also provides loop and iteration functionality which automatically set the
    region Info prefix for each of the solvers returned by the '[]' operator or
    iterator.

Usage
    Given the \c regionSolvers named solvers:
    \verbatim
        // Create the region meshes and solvers
        regionSolvers solvers(runTime);
    \endverbatim

    The list of solvers can be looped over:
    \verbatim
        forAll(solvers, i)
        {
            solvers[i].momentumPredictor();
        }
    \endverbatim
    where the '[]' operator sets the region Info prefix.  After the loop the
    region Info prefix remains set to the last region prefix and so for global
    messages, e.g. the global time-step the Info prefix must be specifically
    reset to spaces by calling the \c setGlobalPrefix() function.

    Alternatively the list of solvers can be iterated over:
    \verbatim
        forAllIter(regionSolvers, solvers, solver)
        {
            solver->momentumPredictor();
        }
    \endverbatim
    where the iterator increment sets the region \c Info prefix and at the end
    automatically resets the \c Info prefix to spaces.

SourceFiles
    regionSolvers.C

\*---------------------------------------------------------------------------*/

#ifndef regionSolvers_H
#define regionSolvers_H

#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionSolvers Declaration
\*---------------------------------------------------------------------------*/

class regionSolvers
{
    // Private Member Data

        //- Region meshes
        PtrList<fvMesh> regions_;

        //- Region solvers
        PtrList<solver> solvers_;

        //- Global space padding prefix
        string prefix0_;

        //- List of space padded region prefixes
        stringList prefixes_;


public:

    // Constructors

        //- Construct from components
        regionSolvers(const Time& runTime);

        //- Disallow default bitwise copy construction
        regionSolvers(const regionSolvers&) = delete;


    //- Destructor
    ~regionSolvers();


    // Member Functions

        //- Return the number of region solvers
        inline label size() const;

        //- Set the Info prefix to space padding for global messages
        void setGlobalPrefix() const;

        //- Set the Info prefix to the space padded region name
        void setPrefix(const label i) const;

        //- Reset the Info prefix to null
        void resetPrefix() const;


    // Iterator

        class iterator;
        friend class iterator;

        // Iterator class to loop over the solvers, setting the prefix for each
        class iterator
        {
            // Private Data

                //- Reference to the regionSolvers
                regionSolvers& regionSolvers_;

                //- Current solver index
                label index_;

        public:

            friend class regionSolvers;

            // Constructors

                //- Construct for the regionSolvers
                inline explicit iterator(regionSolvers&);

            // Member Operators

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline solver& operator*();
                inline solver& operator()();
                inline solver* operator->();

                inline iterator operator++();
                inline iterator operator++(int);
        };

        //- Return an iterator to begin traversing the solvers
        inline iterator begin();

        //- Return an iterator to end traversing the solvers
        inline iterator end();


    // Member Operators

        //- Cast to the list of solvers
        inline operator PtrList<solver>&();

        //- Set the region i prefix and return the corresponding solver
        solver& operator[](const label i);

        //- Disallow default bitwise assignment
        void operator=(const regionSolvers&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolversI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
