/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::SolverPerformance

Description
    SolverPerformance is the class returned by the LduMatrix solver
    containing performance statistics.

SourceFiles
    SolverPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef SolverPerformance_H
#define SolverPerformance_H

#include "word.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class SolverPerformance;

template<class Type>
SolverPerformance<Type> max
(
    const SolverPerformance<Type>&,
    const SolverPerformance<Type>&
);

template<class Type>
Istream& operator>>
(
    Istream&,
    SolverPerformance<Type>&
);

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const SolverPerformance<Type>&
);


/*---------------------------------------------------------------------------*\
                       Class SolverPerformance Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SolverPerformance
{
    word   solverName_;
    word   fieldName_;
    Type   initialResidual_;
    Type   finalResidual_;
    label  noIterations_;
    bool   converged_;
    FixedList<bool, pTraits<Type>::nComponents> singular_;


public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("SolverPerformance");

        //- Large Type for the use in solvers
        static const scalar great_;

        //- Small Type for the use in solvers
        static const scalar small_;

        //- Very small Type for the use in solvers
        static const scalar vsmall_;


    // Constructors

        SolverPerformance()
        :
            initialResidual_(pTraits<Type>::zero),
            finalResidual_(pTraits<Type>::zero),
            noIterations_(0),
            converged_(false),
            singular_(false)
        {}


        SolverPerformance
        (
            const word&  solverName,
            const word&  fieldName,
            const Type&  iRes = pTraits<Type>::zero,
            const Type&  fRes = pTraits<Type>::zero,
            const label  nIter = 0,
            const bool   converged = false,
            const bool   singular = false
        )
        :
            solverName_(solverName),
            fieldName_(fieldName),
            initialResidual_(iRes),
            finalResidual_(fRes),
            noIterations_(nIter),
            converged_(converged),
            singular_(singular)
        {}


    // Member functions

        //- Return solver name
        const word& solverName() const
        {
            return solverName_;
        }

        //- Return solver name
        word& solverName()
        {
            return solverName_;
        }


        //- Return field name
        const word& fieldName() const
        {
            return fieldName_;
        }


        //- Return initial residual
        const Type& initialResidual() const
        {
            return initialResidual_;
        }

        //- Return initial residual
        Type& initialResidual()
        {
            return initialResidual_;
        }


        //- Return final residual
        const Type& finalResidual() const
        {
            return finalResidual_;
        }

        //- Return final residual
        Type& finalResidual()
        {
            return finalResidual_;
        }


        //- Return number of iterations
        label nIterations() const
        {
            return noIterations_;
        }

        //- Return number of iterations
        label& nIterations()
        {
            return noIterations_;
        }


        //- Has the solver converged?
        bool converged() const
        {
            return converged_;
        }

        //- Is the matrix singular?
        bool singular() const;

        //- Check, store and return convergence
        bool checkConvergence
        (
            const Type& tolerance,
            const Type& relTolerance
        );

        //- Singularity test
        bool checkSingularity(const Type& residual);

        //- Print summary of solver performance to the given stream
        void print(Ostream& os) const;


    // Member Operators

        bool operator!=(const SolverPerformance<Type>&) const;


    // Friend functions

        //- Return the element-wise maximum of two SolverPerformance<Type>s
        friend SolverPerformance<Type> max <Type>
        (
            const SolverPerformance<Type>&,
            const SolverPerformance<Type>&
        );


    // Ostream Operator

        friend Istream& operator>> <Type>
        (
            Istream&,
            SolverPerformance<Type>&
        );

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const SolverPerformance<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolverPerformance(Type)                                           \
                                                                              \
typedef Foam::SolverPerformance<Type>                                         \
    solverPerformance##Type;                                                  \
                                                                              \
defineNamedTemplateTypeNameAndDebug(solverPerformance##Type, 0);              \
                                                                              \
template<>                                                                    \
const scalar solverPerformance##Type::great_(1e20);                           \
                                                                              \
template<>                                                                    \
const scalar solverPerformance##Type::small_(1e-20);                          \
                                                                              \
template<>                                                                    \
const scalar solverPerformance##Type::vsmall_(VSMALL);                        \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SolverPerformance.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
