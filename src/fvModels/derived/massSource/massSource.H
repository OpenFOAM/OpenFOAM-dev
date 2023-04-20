/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

    If the mass flow rate is positive then user-supplied fixed property values
    are introduced to the field equations. If the mass flow rate is negative
    then properties are removed at their current value.

Usage
    Example usage:
    \verbatim
    massSource
    {
        type            massSource;

        select          cellSet;
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

    If the mass flow rate is positive then values should be provided for all
    solved for fields. Warnings will be issued if values are not provided for
    fields for which transport equations are solved. Warnings will also be
    issued if values are provided for fields which are not solved for.

SourceFiles
    massSource.C

\*---------------------------------------------------------------------------*/

#ifndef massSource_H
#define massSource_H

#include "fvModel.H"
#include "fvCellSet.H"
#include "HashPtrTable.H"
#include "unknownTypeFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class massSourceBase Declaration
\*---------------------------------------------------------------------------*/

class massSourceBase
:
    public fvModel
{
private:

    // Private Data

        //- Name of the phase
        word phaseName_;

        //- Name of the density field
        word rhoName_;

        //- Name of the energy field
        word heName_;

        //- Name of the temperature field
        word TName_;

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Field values
        HashPtrTable<unknownTypeFunction1> fieldValues_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the mass flow rate
        virtual scalar massFlowRate() const = 0;


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


protected:

    // Protected Member Functions

        //- Read the set
        void readSet();

        //- Read the field values
        void readFieldValues();


public:

    //- Runtime type information
    TypeName("massSourceBase");


    // Constructors

        //- Construct from explicit source name and mesh
        massSourceBase
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        massSourceBase(const massSourceBase&) = delete;


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

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const massSourceBase&) = delete;
};


/*---------------------------------------------------------------------------*\
                          Class massSource Declaration
\*---------------------------------------------------------------------------*/

class massSource
:
    public massSourceBase
{
private:

    // Private Data

        //- Mass flow rate
        autoPtr<Function1<scalar>> massFlowRate_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the mass flow rate
        virtual scalar massFlowRate() const;


public:

    //- Runtime type information
    TypeName("massSource");


    // Constructors

        //- Construct from explicit source name and mesh
        massSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
