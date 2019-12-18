/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::ODESolver

Description
    Abstract base-class for ODE system solvers

SourceFiles
    ODESolver.C

\*---------------------------------------------------------------------------*/

#ifndef ODESolver_H
#define ODESolver_H

#include "ODESystem.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ODESolver Declaration
\*---------------------------------------------------------------------------*/

class ODESolver
{

protected:

    // Protected data

        //- Reference to ODESystem
        const ODESystem& odes_;

        //- Maximum size of the ODESystem
        const label maxN_;

        //- Size of the ODESystem (adjustable)
        mutable label n_;

        //- Absolute convergence tolerance per step
        scalarField absTol_;

        //- Relative convergence tolerance per step
        scalarField relTol_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;


    // Protected Member Functions

        //- Return the normalized scalar error
        scalar normalizeError
        (
            const scalarField& y0,
            const scalarField& y,
            const scalarField& err
        ) const;


public:

    friend class ODESystem;

    //- Runtime type information
    TypeName("ODESolver");

    class stepState
    {
        public:

        const bool forward;
        scalar dxTry;
        scalar dxDid;
        bool first;
        bool last;
        bool reject;
        bool prevReject;

        stepState(const scalar dx)
        :
            forward(dx > 0 ? true : false),
            dxTry(dx),
            dxDid(0),
            first(true),
            last(false),
            reject(false),
            prevReject(false)
        {}
    };


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ODESolver,
            dictionary,
            (const ODESystem& ode, const dictionary& dict),
            (ode, dict)
        );


    // Constructors

        //- Construct for given ODESystem
        ODESolver(const ODESystem& ode, const dictionary& dict);

        //- Construct for given ODESystem specifying tolerances
        ODESolver
        (
            const ODESystem& ode,
            const scalarField& absTol,
            const scalarField& relTol
        );

        //- Disallow default bitwise copy construction
        ODESolver(const ODESolver&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<ODESolver> New
        (
            const ODESystem& ode,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ODESolver()
    {}


    // Member Functions

        //- Return the number of equations to solve
        inline label nEqns() const;

        //- Return access to the absolute tolerance field
        inline scalarField& absTol();

        //- Return access to the relative tolerance field
        inline scalarField& relTol();

        //- Resize the ODE solver
        virtual bool resize() = 0;

        template<class Type>
        static inline void resizeField(UList<Type>& f, const label n);

        template<class Type>
        inline void resizeField(UList<Type>& f) const;

        inline void resizeMatrix(scalarSquareMatrix& m) const;

        //- Solve the ODE system from the current state xStart, y
        //  and the optional index into the list of systems to solve li
        //  as far as possible up to dxTry adjusting the step as necessary
        //  to provide a solution within the specified tolerance.
        //  Update the state and return an estimate for the next step in dxTry
        virtual void solve
        (
            scalar& x,
            scalarField& y,
            const label li,
            scalar& dxTry
        ) const;

        //- Solve the ODE system from the current state xStart, y
        //  and the optional index into the list of systems to solve li
        //  as far as possible up to step.dxTry adjusting the step as necessary
        //  to provide a solution within the specified tolerance.
        //  Update the state and return the step actually taken in step.dxDid
        //  and an estimate for the next step in step.dxTry
        virtual void solve
        (
            scalar& x,
            scalarField& y,
            const label li,
            stepState& step
        ) const;

        //- Solve the ODE system from the current state xStart, y
        //  and the optional index into the list of systems to solve li
        //  to xEnd and return an estimate for the next step in dxTry
        virtual void solve
        (
            const scalar xStart,
            const scalar xEnd,
            scalarField& y,
            const label li,
            scalar& dxEst
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ODESolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ODESolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
