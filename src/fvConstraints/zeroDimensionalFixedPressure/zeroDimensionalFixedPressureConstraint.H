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
    Foam::fv::zeroDimensionalFixedPressureConstraint

Description
    Zero-dimensional fixed pressure constraint. Should be used in conjunction
    with the zeroDimensionalFixedPressureModel.

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

        // Name of the pressure field, default = p
        //p               p;

        // Name of the density field, default = rho
        //rho             rho;

        // Pressure value
        pressure        1e5;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef zeroDimensionalFixedPressureConstraint_H
#define zeroDimensionalFixedPressureConstraint_H

#include "fvConstraint.H"
#include "Function1.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class zeroDimensionalFixedPressureModel;

/*---------------------------------------------------------------------------*\
             Class zeroDimensionalFixedPressureConstraint Declaration
\*---------------------------------------------------------------------------*/

class zeroDimensionalFixedPressureConstraint
:
    public fvConstraint
{
    // Private data

        //- Pressure field name, default = p
        word pName_;

        //- Density field name, default = rho
        word rhoName_;

        //- The pressure value
        autoPtr<Function1<scalar>> p_;

        //- The mass or volume source
        mutable autoPtr<volScalarField::Internal> sourcePtr_;


    // Private member functions

        //- Access the corresponding model
        const zeroDimensionalFixedPressureModel& model() const;

        //- Get the mass source
        template<class AlphaFieldType>
        tmp<volScalarField::Internal> massSource
        (
            const AlphaFieldType& alpha,
            const volScalarField::Internal& rho
        ) const;

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("zeroDimensionalFixedPressure");


    // Constructors

        //- Construct from dictionary
        zeroDimensionalFixedPressureConstraint
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zeroDimensionalFixedPressureConstraint();


    // Member Functions

        // Access

            //- Pressure field name
            inline const word& pName() const
            {
                return pName_;
            }

            //- Density field name
            inline const word& rhoName() const
            {
                return rhoName_;
            }


        // Checks

            //- Return the list of fields constrained by the fvConstraint
            virtual wordList constrainedFields() const;


        // Constraints

            //- Return the mass or volume source for the pressure equation
            const volScalarField::Internal& pEqnSource
            (
                fvMatrix<scalar>& pEqn
            ) const;

            //- Return the mass source
            tmp<volScalarField::Internal> massSource
            (
                const volScalarField::Internal& rho
            ) const;

            //- Return the mass source for a given phase
            tmp<volScalarField::Internal> massSource
            (
                const volScalarField::Internal& alpha,
                const volScalarField::Internal& rho
            ) const;

            //- Apply the constraint to the pressure equation
            virtual bool constrain
            (
                fvMatrix<scalar>& pEqn,
                const word& fieldName
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
