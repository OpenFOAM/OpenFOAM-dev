/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::fv::massTransfer

Description
    Base class for mass transfers between phases

SourceFiles
    massTransfer.C
    massTransferI.H

\*---------------------------------------------------------------------------*/

#ifndef massTransfer_H
#define massTransfer_H

#include "fvSpecificSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class massTransfer Declaration
\*---------------------------------------------------------------------------*/

class massTransfer
:
    public fvSpecificSource
{
private:

    // Private Data

        //- Names of the phases
        const Pair<word> phaseNames_;

        //- Names of the volume fraction fields
        const Pair<word> alphaNames_;

        //- Names of the density fields
        const Pair<word> rhoNames_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


protected:

    // Protected Member Functions

        //- Lookup the phase names
        const Pair<word> lookupPhaseNames(const dictionary& dict) const;

        //- Lookup the phase field names
        const Pair<word> lookupPhaseFieldNames
        (
            const dictionary& dict,
            const word& name
        ) const;

        //- Get the index of the name in the given pair, or -1
        static inline label index(const Pair<word>& names, const word& name);

        //- Get the sign associated with the name in the given pair; i.e., -1
        //  if it is the first name, +1 if the second, and zero if neither
        static inline label sign(const Pair<word>& names, const word& name);

        //- Return the density
        tmp<volScalarField::Internal> rho(const label i) const;

        //- Return the corresponding field in the other phase, or a null
        //  reference if no such field exists
        template<class Type>
        const VolField<Type>& otherField(const VolField<Type>& field) const;


        // Sources

            //- Add a source term to an equation
            template<class Type>
            void addSupType
            (
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a scalar equation
            void addSupType
            (
                const volScalarField& alphaOrField,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& alphaOrRho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a compressible scalar equation
            void addSupType
            (
                const volScalarField& alphaOrRho,
                const volScalarField& rhoOrField,
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
    TypeName("massTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        massTransfer
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Access

            //- Return the names of the phases
            inline const Pair<word>& phaseNames() const;

            //- Return the names of the volume fraction fields
            inline const Pair<word>& alphaNames() const;

            //- Return the names of the density fields
            inline const Pair<word>& rhoNames() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const = 0;

            //- Return the source value
            virtual tmp<DimensionedField<scalar, volMesh>> S
            (
                const word& fieldName
            ) const;

            //- Add a source term to a field-less proxy equation
            virtual void addSup(fvMatrix<scalar>& eqn) const;

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_FIELD_SUP);

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_FIELD_SUP);

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP);


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "massTransferI.H"

#ifdef NoRepository
    #include "massTransferTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
