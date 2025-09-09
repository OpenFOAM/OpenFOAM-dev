/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::Lagrangian::constantFluxCarrierTransfer

Description
    Model to represent mass transfer into the carrier phase with a constant
    surface volume or mass flux. This model provides a crude approximation of
    phase change which neglects any thermal effects. It can be used by clouds
    and in solvers that do not feature any thermodynamic modelling.

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        constantFluxCarrierTransfer;
        volumeFlux  1 [ml/cm^2/s];
    }
    \endverbatim

SourceFiles
    constantFluxCarrierTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef constantFluxCarrierTransfer_H
#define constantFluxCarrierTransfer_H

#include "cloudLagrangianModel.H"
#include "LagrangianSource.H"
#include "sharedRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                    Class constantFluxCarrierTransfer Declaration
\*---------------------------------------------------------------------------*/

class constantFluxCarrierTransfer
:
    public LagrangianSource,
    private cloudLagrangianModel
{
    // Private Member Data

        //- Volume or mass flux
        dimensionedScalar flux_;

        //- Sum of the implicit coefficients of the mass transfer models
        sharedRegIOobject<LagrangianSubScalarField> sumDeltaTSp_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& modelDict);

        //- Return the implicit coefficient for the transfer
        tmp<LagrangianSubScalarField> Sp(const LagrangianSubMesh&) const;

        //- Add a source term to the volume or mass equation
        void addSupType
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& vOrM,
            LagrangianEqn<scalar>& eqn
        ) const;

        //- Add a source term to an equation
        template<class Type>
        void addSupType
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubSubField<Type>& field,
            LagrangianEqn<Type>& eqn
        ) const;

        //- Add a source term to a volume or mass-weighted equation
        template<class Type>
        void addSupType
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& vOrM,
            const LagrangianSubSubField<Type>& field,
            LagrangianEqn<Type>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("constantFluxCarrierTransfer");


    // Constructors

        //- Construct from components
        constantFluxCarrierTransfer
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Return the name of the volume or mass field
        virtual wordList addSupFields() const;

        //- Return true if this is a Lagrangian field or a field of the
        //  Eulerian carrier phase
        virtual bool addsSupToField
        (
            const word& fieldName,
            const word& eqnFieldName
        ) const;

        //- Update the transfer rate and remove any consumed particles
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Return the source value
        virtual tmp<LagrangianSubScalarField> source
        (
            const word& fieldName,
            const LagrangianSubMesh& subMesh
        ) const;

        //- Hook before source evaluation
        virtual void preAddSup
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Add a source term to an equation
        FOR_ALL_FIELD_TYPES(DEFINE_LAGRANGIAN_MODEL_ADD_FIELD_SUP)

        //- Add a source term to a volume or mass-weighted equation
        FOR_ALL_FIELD_TYPES(DEFINE_LAGRANGIAN_MODEL_ADD_V_OR_M_FIELD_SUP)

        //- Hook after source evaluation
        virtual void postAddSup
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Read dictionary
        virtual bool read(const dictionary& modelDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
