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
    Foam::LduMatrix

Description
    LduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.

    Addressing arrays must be supplied for the upper and lower triangles.

Note
    It might be better if this class were organised as a hierarchy starting
    from an empty matrix, then deriving diagonal, symmetric and asymmetric
    matrices.

SourceFiles
    LduMatrixATmul.C
    LduMatrix.C
    LduMatrixOperations.C
    LduMatrixSolver.C
    LduMatrixPreconditioner.C
    LduMatrixTests.C
    LduMatrixUpdateMatrixInterfaces.C

\*---------------------------------------------------------------------------*/

#ifndef LduMatrix_H
#define LduMatrix_H

#include "lduMesh.H"
#include "Field.H"
#include "FieldField.H"
#include "LduInterfaceFieldPtrsList.H"
#include "SolverPerformance.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type, class DType, class LUType>
class LduMatrix;

template<class Type, class DType, class LUType>
Ostream& operator<<
(
    Ostream&,
    const LduMatrix<Type, DType, LUType>&
);


/*---------------------------------------------------------------------------*\
                          Class LduMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class LduMatrix
{
    // private data

        //- LDU mesh reference
        const lduMesh& lduMesh_;

        //- Diagonal coefficients
        Field<DType> *diagPtr_;

        //- Off-diagonal coefficients
        Field<LUType> *upperPtr_, *lowerPtr_;

        //- Source
        Field<Type> *sourcePtr_;

        //- Field interfaces (processor patches etc.)
        LduInterfaceFieldPtrsList<Type> interfaces_;

        //- Off-diagonal coefficients for interfaces
        FieldField<Field, LUType> interfacesUpper_, interfacesLower_;


public:

    friend class SolverPerformance<Type>;

    //- Abstract base-class for LduMatrix solvers
    class solver
    {
    protected:

        // Protected data

            word fieldName_;
            const LduMatrix<Type, DType, LUType>& matrix_;

            //- Dictionary of controls
            dictionary controlDict_;

            //- Default maximum number of iterations in the solver
            static const label defaultMaxIter_ = 1000;

            //- Maximum number of iterations in the solver
            label maxIter_;

            //- Minimum number of iterations in the solver
            label minIter_;

            //- Final convergence tolerance
            Type tolerance_;

            //- Convergence tolerance relative to the initial
            Type relTol_;


        // Protected Member Functions

            //- Read a control parameter from controlDict
            template<class T>
            inline void readControl
            (
                const dictionary& controlDict,
                T& control,
                const word& controlName
            );


            //- Read the control parameters from the controlDict_
            virtual void readControls();


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                symMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix,
                    const dictionary& solverDict
                ),
                (
                    fieldName,
                    matrix,
                    solverDict
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                asymMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix,
                    const dictionary& solverDict
                ),
                (
                    fieldName,
                    matrix,
                    solverDict
                )
            );


        // Constructors

            solver
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix,
                const dictionary& solverDict
            );


        // Selectors

            //- Return a new solver
            static autoPtr<solver> New
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix,
                const dictionary& solverDict
            );


        // Destructor

            virtual ~solver()
            {}


        // Member Functions

            // Access

                const word& fieldName() const
                {
                    return fieldName_;
                }

                const LduMatrix<Type, DType, LUType>& matrix() const
                {
                    return matrix_;
                }


            //- Read and reset the solver parameters from the given dictionary
            virtual void read(const dictionary& solverDict);

            virtual SolverPerformance<Type> solve
            (
                Field<Type>& psi
            ) const = 0;

            //- Return the matrix norm used to normalise the residual for the
            //  stopping criterion
            Type normFactor
            (
                const Field<Type>& psi,
                const Field<Type>& Apsi,
                Field<Type>& tmpField
            ) const;
    };


    //- Abstract base-class for LduMatrix smoothers
    class smoother
    {
    protected:

        // Protected data

            word fieldName_;
            const LduMatrix<Type, DType, LUType>& matrix_;


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                symMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix
                ),
                (
                    fieldName,
                    matrix
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                asymMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix
                ),
                (
                    fieldName,
                    matrix
                )
            );


        // Constructors

            smoother
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix
            );


        // Selectors

            //- Return a new smoother
            static autoPtr<smoother> New
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix,
                const dictionary& smootherDict
            );


        // Destructor

            virtual ~smoother()
            {}


        // Member Functions

            // Access

                const word& fieldName() const
                {
                    return fieldName_;
                }

                const LduMatrix<Type, DType, LUType>& matrix() const
                {
                    return matrix_;
                }


            //- Smooth the solution for a given number of sweeps
            virtual void smooth
            (
                Field<Type>& psi,
                const label nSweeps
            ) const = 0;
    };


    //- Abstract base-class for LduMatrix preconditioners
    class preconditioner
    {
    protected:

        // Protected data

            //- Reference to the base-solver this preconditioner is used with
            const solver& solver_;


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                symMatrix,
                (
                    const solver& sol,
                    const dictionary& preconditionerDict
                ),
                (sol, preconditionerDict)
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                asymMatrix,
                (
                    const solver& sol,
                    const dictionary& preconditionerDict
                ),
                (sol, preconditionerDict)
            );


        // Constructors

            preconditioner
            (
                const solver& sol
            )
            :
                solver_(sol)
            {}


        // Selectors

            //- Return a new preconditioner
            static autoPtr<preconditioner> New
            (
                const solver& sol,
                const dictionary& preconditionerDict
            );


        // Destructor

            virtual ~preconditioner()
            {}


        // Member Functions

            //- Read and reset the preconditioner parameters
            //  from the given dictionary
            virtual void read(const dictionary& preconditionerDict)
            {}

            //- Return wA the preconditioned form of residual rA
            virtual void precondition
            (
                Field<Type>& wA,
                const Field<Type>& rA
            ) const = 0;

            //- Return wT the transpose-matrix preconditioned form of
            //  residual rT.
            //  This is only required for preconditioning asymmetric matrices.
            virtual void preconditionT
            (
                Field<Type>& wT,
                const Field<Type>& rT
            ) const
            {
                NotImplemented;
            }
    };


    // Static data

        // Declare name of the class and its debug switch
        ClassName("LduMatrix");


    // Constructors

        //- Construct given an LDU addressed mesh.
        //  The coefficients are initially empty for subsequent setting.
        LduMatrix(const lduMesh&);

        //- Copy constructor
        LduMatrix(const LduMatrix<Type, DType, LUType>&);

        //- Copy constructor or re-use as specified.
        LduMatrix(LduMatrix<Type, DType, LUType>&, bool reuse);

        //- Construct given an LDU addressed mesh and an Istream
        //  from which the coefficients are read
        LduMatrix(const lduMesh&, Istream&);


    // Destructor

        ~LduMatrix();


    // Member Functions

        // Access to addressing

            //- Return the LDU mesh from which the addressing is obtained
            const lduMesh& mesh() const
            {
                return lduMesh_;
            }

            //- Return the LDU addressing
            const lduAddressing& lduAddr() const
            {
                return lduMesh_.lduAddr();
            }

            //- Return the patch evaluation schedule
            const lduSchedule& patchSchedule() const
            {
                return lduAddr().patchSchedule();
            }

            //- Return interfaces
            const LduInterfaceFieldPtrsList<Type>& interfaces() const
            {
                return interfaces_;
            }

            //- Return interfaces
            LduInterfaceFieldPtrsList<Type>& interfaces()
            {
                return interfaces_;
            }


        // Access to coefficients

            Field<DType>& diag();
            Field<LUType>& upper();
            Field<LUType>& lower();
            Field<Type>& source();

            FieldField<Field, LUType>& interfacesUpper()
            {
                return interfacesUpper_;
            }

            FieldField<Field, LUType>& interfacesLower()
            {
                return interfacesLower_;
            }


            const Field<DType>& diag() const;
            const Field<LUType>& upper() const;
            const Field<LUType>& lower() const;
            const Field<Type>& source() const;

            const FieldField<Field, LUType>& interfacesUpper() const
            {
                return interfacesUpper_;
            }

            const FieldField<Field, LUType>& interfacesLower() const
            {
                return interfacesLower_;
            }


            bool hasDiag() const
            {
                return (diagPtr_);
            }

            bool hasUpper() const
            {
                return (upperPtr_);
            }

            bool hasLower() const
            {
                return (lowerPtr_);
            }

            bool hasSource() const
            {
                return (sourcePtr_);
            }

            bool diagonal() const
            {
                return
                (
                    diagPtr_
                 && Pstream::parRun()
                    ?
                        !lowerPtr_ && !upperPtr_
                    :
                        !(lowerPtr_ && lowerPtr_->size())
                     && !(upperPtr_ && upperPtr_->size())
                );
            }

            bool symmetric() const
            {
                return
                (
                    diagPtr_
                 && Pstream::parRun()
                    ?
                        !lowerPtr_ && upperPtr_
                    :
                        !(lowerPtr_ && lowerPtr_->size())
                     && (upperPtr_ && upperPtr_->size())
                );
            }

            bool asymmetric() const
            {
                return
                (
                    diagPtr_
                 && Pstream::parRun()
                    ?
                        lowerPtr_ && upperPtr_
                    :
                        (lowerPtr_ && lowerPtr_->size())
                     && (upperPtr_ && upperPtr_->size())
                );
            }


        // operations

            void sumDiag();
            void negSumDiag();

            void sumMagOffDiag(Field<LUType>& sumOff) const;

            //- Matrix multiplication
            void Amul(Field<Type>&, const tmp<Field<Type>>&) const;

            //- Matrix transpose multiplication
            void Tmul(Field<Type>&, const tmp<Field<Type>>&) const;


            //- Sum the coefficients on each row of the matrix
            void sumA(Field<Type>&) const;


            void residual(Field<Type>& rA, const Field<Type>& psi) const;

            tmp<Field<Type>> residual(const Field<Type>& psi) const;


            //- Initialise the update of interfaced interfaces
            //  for matrix operations
            void initMatrixInterfaces
            (
                const FieldField<Field, LUType>& interfaceCoeffs,
                const Field<Type>& psiif,
                Field<Type>& result
            ) const;

            //- Update interfaced interfaces for matrix operations
            void updateMatrixInterfaces
            (
                const FieldField<Field, LUType>& interfaceCoeffs,
                const Field<Type>& psiif,
                Field<Type>& result
            ) const;


            tmp<Field<Type>> H(const Field<Type>&) const;
            tmp<Field<Type>> H(const tmp<Field<Type>>&) const;

            tmp<Field<Type>> faceH(const Field<Type>&) const;
            tmp<Field<Type>> faceH(const tmp<Field<Type>>&) const;


    // Member Operators

        void operator=(const LduMatrix<Type, DType, LUType>&);

        void negate();

        void operator+=(const LduMatrix<Type, DType, LUType>&);
        void operator-=(const LduMatrix<Type, DType, LUType>&);

        void operator*=(const scalarField&);
        void operator*=(scalar);


    // Ostream operator

        friend Ostream& operator<< <Type, DType, LUType>
        (
            Ostream&,
            const LduMatrix<Type, DType, LUType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLduMatrix(Type, DType, LUType)                                     \
                                                                               \
typedef Foam::LduMatrix<Type, DType, LUType>                                   \
    ldu##Type##DType##LUType##Matrix;                                          \
                                                                               \
defineNamedTemplateTypeNameAndDebug(ldu##Type##DType##LUType##Matrix, 0);      \
                                                                               \
                                                                               \
typedef LduMatrix<Type, DType, LUType>::smoother                               \
    ldu##Type##DType##LUType##Smoother;                                        \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Smoother,                                        \
    symMatrix                                                                  \
);                                                                             \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Smoother,                                        \
    asymMatrix                                                                 \
);                                                                             \
                                                                               \
                                                                               \
typedef LduMatrix<Type, DType, LUType>::preconditioner                         \
    ldu##Type##DType##LUType##Preconditioner;                                  \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Preconditioner,                                  \
    symMatrix                                                                  \
);                                                                             \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Preconditioner,                                  \
    asymMatrix                                                                 \
);                                                                             \
                                                                               \
                                                                               \
typedef LduMatrix<Type, DType, LUType>::solver                                 \
    ldu##Type##DType##LUType##Solver;                                          \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Solver,                                          \
    symMatrix                                                                  \
);                                                                             \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Solver,                                          \
    asymMatrix                                                                 \
);


#define makeLduPreconditioner(Precon, Type, DType, LUType)                     \
                                                                               \
typedef Precon<Type, DType, LUType>                                            \
    Precon##Type##DType##LUType##Preconditioner;                               \
defineNamedTemplateTypeNameAndDebug                                            \
(                                                                              \
    Precon##Type##DType##LUType##Preconditioner,                               \
    0                                                                          \
);

#define makeLduSymPreconditioner(Precon, Type, DType, LUType)                  \
                                                                               \
LduMatrix<Type, DType, LUType>::preconditioner::                               \
addsymMatrixConstructorToTable<Precon##Type##DType##LUType##Preconditioner>    \
add##Precon##Type##DType##LUType##PreconditionerSymMatrixConstructorToTable_;

#define makeLduAsymPreconditioner(Precon, Type, DType, LUType)                 \
                                                                               \
LduMatrix<Type, DType, LUType>::preconditioner::                               \
addasymMatrixConstructorToTable<Precon##Type##DType##LUType##Preconditioner>   \
add##Precon##Type##DType##LUType##PreconditionerAsymMatrixConstructorToTable_;


#define makeLduSmoother(Smoother, Type, DType, LUType)                         \
                                                                               \
typedef Smoother<Type, DType, LUType>                                          \
    Smoother##Type##DType##LUType##Smoother;                                   \
                                                                               \
defineNamedTemplateTypeNameAndDebug                                            \
(                                                                              \
    Smoother##Type##DType##LUType##Smoother,                                   \
    0                                                                          \
);

#define makeLduSymSmoother(Smoother, Type, DType, LUType)                      \
                                                                               \
LduMatrix<Type, DType, LUType>::smoother::                                     \
    addsymMatrixConstructorToTable<Smoother##Type##DType##LUType##Smoother>    \
    add##Smoother##Type##DType##LUType##SymMatrixConstructorToTable_;

#define makeLduAsymSmoother(Smoother, Type, DType, LUType)                     \
                                                                               \
LduMatrix<Type, DType, LUType>::smoother::                                     \
    addasymMatrixConstructorToTable<Smoother##Type##DType##LUType##Smoother>   \
    add##Smoother##Type##DType##LUType##AsymMatrixConstructorToTable_;


#define makeLduSolver(Solver, Type, DType, LUType)                             \
                                                                               \
typedef Solver<Type, DType, LUType>                                            \
    Solver##Type##DType##LUType##Solver;                                       \
                                                                               \
defineNamedTemplateTypeNameAndDebug                                            \
(                                                                              \
    Solver##Type##DType##LUType##Solver,                                       \
    0                                                                          \
);

#define makeLduSymSolver(Solver, Type, DType, LUType)                          \
                                                                               \
LduMatrix<Type, DType, LUType>::solver::                                       \
    addsymMatrixConstructorToTable<Solver##Type##DType##LUType##Solver>        \
    add##Solver##Type##DType##LUType##SymMatrixConstructorToTable_;

#define makeLduAsymSolver(Solver, Type, DType, LUType)                         \
                                                                               \
LduMatrix<Type, DType, LUType>::solver::                                       \
    addasymMatrixConstructorToTable<Solver##Type##DType##LUType##Solver>       \
    add##Solver##Type##DType##LUType##AsymMatrixConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LduMatrixI.H"
    #include "LduMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
