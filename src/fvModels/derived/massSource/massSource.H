/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fv::massSource

Description
    This fvModel applies a mass source to the continuity equation and to all
    field equations.

Usage
    Example usage:
    \verbatim
    massSource
    {
        type            massSource;

        selectionMode   cellSet;
        cellSet         massSource;

        massFlowRate    1e-4;

        fieldValues
        {
            U               (10 0 0);
            T               300;
            k               0.375;
            epsilon         14.855;
        }
    }
    \endverbatim

    Values should be provided for all solved for fields. Warnings will be
    issued if values are not provided for fields for which transport equations
    are solved. Warnings will also be issued if values are provided for fields
    which are not solved for.

SourceFiles
    massSource.C

\*---------------------------------------------------------------------------*/

#ifndef massSource_H
#define massSource_H

#include "fvModel.H"
#include "fvCellSet.H"
#include "HashPtrTable.H"
#include "objectFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class massSource Declaration
\*---------------------------------------------------------------------------*/

class massSource
:
    public fvModel
{
    // Private Data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Name of the phase
        word phaseName_;

        //- Name of the density field
        word rhoName_;

        //- Name of the energy field
        word heName_;

        //- Name of the temperature field
        word TName_;

        //- Mass source
        autoPtr<Function1<scalar>> massFlowRate_;

        //- Field values
        HashPtrTable<objectFunction1> fieldValues_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Read field values for a given primitive type
        template<class Type>
        void readFieldValues();


        // Sources

            //- Add a source term to an equation
            template<class Type>
            void addGeneralSupType
            (
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;

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
    TypeName("massSource");


    // Constructors

        //- Construct from explicit source name and mesh
        massSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        massSource(const massSource&) = delete;


    // Member Functions

        // Access

            //- Return the name of the density field
            const word& rhoName() const
            {
                return rhoName_;
            }


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


        // Mesh motion

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const massSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
