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
    Foam::Lagrangian::collisionPhaseTransfer

Description
    Model to represent the absorption of droplets or bubbles into a
    corresponding Eulerian phase as a result of collisions with that phase (and
    potentially other non-carrier Eulerian phases). This can be used to
    represent (e.g.) a spray of droplets hitting and being absorbed by a VoF
    interface. The collision probability is calculated as proportional to the
    Eulerian phase's surface area per unit volume and the distance travelled by
    the particle.

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        collisionPhaseTransfer;
    }
    \endverbatim

SourceFiles
    collisionPhaseTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef collisionPhaseTransfer_H
#define collisionPhaseTransfer_H

#include "cloudLagrangianModel.H"
#include "LagrangianSource.H"
#include "CarrierField.H"
#include "restartableRandomGenerator.H"
#include "sharedRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                    Class collisionPhaseTransfer Declaration
\*---------------------------------------------------------------------------*/


class collisionPhaseTransfer
:
    public LagrangianSource,
    private cloudLagrangianModel
{
    // Private Member Data

        //- Carrier phase fraction
        const CarrierField<scalar>& alphac_;

        //- Random generator
        restartableRandomGenerator rndGen_;

        //- Cached hit probability for this step
        autoPtr<LagrangianSubScalarField> PhitPtr_;

        //- Sum of the implicit coefficients of the mass transfer models
        sharedRegIOobject<LagrangianSubScalarField> sumDeltaTSp_;


    // Private Member Functions

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
    TypeName("collisionPhaseTransfer");


    // Constructors

        //- Construct from components
        collisionPhaseTransfer
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
        //  corresponding Eulerian phase
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

        //- Add a fractional source term
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            LagrangianEqn<scalar>& eqn
        ) const;

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

        //- Write state
        virtual void writeProcessorState(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
