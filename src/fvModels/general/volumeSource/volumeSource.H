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
    Foam::fv::volumeSource

Description
    This fvModel applies a volume source to the continuity equation and to all
    field equations. It can be applied to incompressible solvers, such as
    incompressibleFluid and incompressibleVoF. For compressible solvers, use
    the massSource model instead.

    If the volumetric flow rate is positive then user-supplied fixed property
    values are introduced to the field equations. If the volumetric flow rate
    is negative then properties are removed at their current value.

Usage
    Example usage:
    \verbatim
    volumeSource
    {
        type            volumeSource;

        select          cellSet;
        cellSet         volumeSource;

        volumetricFlowRate 1e-4;

        fieldValues
        {
            U               (10 0 0);
            k               0.375;
            epsilon         14.855;
        }
    }
    \endverbatim

    If the volumetric flow rate is positive then values should be provided for
    all solved for fields. Warnings will be issued if values are not provided
    for fields for which transport equations are solved. Warnings will also be
    issued if values are provided for fields which are not solved for.

SourceFiles
    volumeSource.C

See also
    Foam::fv::massSource

\*---------------------------------------------------------------------------*/

#ifndef volumeSource_H
#define volumeSource_H

#include "fvTotalSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class volumeSource Declaration
\*---------------------------------------------------------------------------*/

class volumeSource
:
    public fvTotalSource
{
private:

    // Private Data

        //- Name of the volume fraction field (if any)
        word alphaName_;

        //- The set of cells the source applies to
        autoPtr<fvCellSet> setPtr_;

        //- Volumetric flow rate
        autoPtr<Function1<scalar>> volumetricFlowRate_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


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
    TypeName("volumeSource");


    // Constructors

        //- Construct from explicit source name and mesh
        volumeSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return the cells that the source applies to
            virtual labelUList cells() const;

            //- Return the number of cells that the source applies to
            virtual label nCells() const;

            //- Return the volume of cells that the source applies to
            virtual scalar V() const;


        // Sources

            //- Return the source value
            virtual dimensionedScalar S() const;

            //- Add a source term to a field-less proxy equation
            virtual void addSup(fvMatrix<scalar>& eqn) const;

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
