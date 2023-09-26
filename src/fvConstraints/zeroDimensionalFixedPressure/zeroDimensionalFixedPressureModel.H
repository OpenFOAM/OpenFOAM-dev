/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::fv::zeroDimensionalFixedPressureModel

Description
    Zero-dimensional fixed pressure source. Should be used in conjunction
    with the zeroDimensionalFixedPressureConstraint.

    This constraint and model facilitates specification of a constant or
    time-varying pressure. It adds mass source terms proportional to the error
    that remains when the pressure equation is evaluated at the desired
    pressure. Iteration may be necessary to converge the constraint in the case
    of non-linear equations of state.

    Properties are added or removed with their current value. The model
    therefore represents a uniform expansion or contraction in infinite space.

Usage
    Example usage:
    \verbatim
    {
        type            zeroDimensionalFixedPressure;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef zeroDimensionalFixedPressureModel_H
#define zeroDimensionalFixedPressureModel_H

#include "fvModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class zeroDimensionalFixedPressureConstraint;

/*---------------------------------------------------------------------------*\
              Class zeroDimensionalFixedPressureModel Declaration
\*---------------------------------------------------------------------------*/

class zeroDimensionalFixedPressureModel
:
    public fvModel
{
    // Private member functions

        //- Access the corresponding constraint
        const zeroDimensionalFixedPressureConstraint& constraint() const;


        // Sources

            //- Add a source term to an equation
            template<class Type>
            void addSupType
            (
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a compressible continuity equation
            void addSupType
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a phase continuity equation
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add a source term to a phase equation
            template<class Type>
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;


public:

    //- Runtime type information
    TypeName("zeroDimensionalFixedPressure");


    // Constructors

        //- Construct from dictionary
        zeroDimensionalFixedPressureModel
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zeroDimensionalFixedPressureModel();


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Sources

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_FIELD_SUP)

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_FIELD_SUP)

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP)


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
