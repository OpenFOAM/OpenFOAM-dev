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
    Foam::fv::zeroDimensional::constraintSource

Description
    Base class for zero-dimensional constraint sources.

\*---------------------------------------------------------------------------*/

#ifndef constraintSource_H
#define constraintSource_H

#include "zeroDimensionalFvModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace zeroDimensional
{

/*---------------------------------------------------------------------------*\
                      Class constraintSource Declaration
\*---------------------------------------------------------------------------*/

class constraintSource
:
    public zeroDimensionalFvModel
{
    // Private member functions

        // Sources

            //- Return the mass source
            virtual tmp<volScalarField::Internal> dmdt() const = 0;

            //- Add a source term to an equation
            template<class Type>
            void addSupType(fvMatrix<Type>& eqn, const word& fieldName) const;

            //- Add a source term to a scalar equation
            void addSupType(fvMatrix<scalar>& eqn, const word& fieldName) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;

            //- Add a source term to a scalar compressible equation
            void addSupType
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add a source term to a phase equation
            template<class Type>
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;


public:

    //- Runtime type information
    TypeName("constraintSource");


    // Constructors

        //- Construct from dictionary
        constraintSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constraintSource();


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_SUP);

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_SUP);

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_SUP);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zeroDimensional
} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
