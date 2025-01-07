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
    Foam::turbulentDispersion

Description
    ...

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        turbulentDispersion;
    }
    \endverbatim

SourceFiles
    turbulentDispersion.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDispersion_H
#define turbulentDispersion_H

#include "cloudLagrangianModel.H"
#include "drag.H"
#include "momentumTransportModel.H"
#include "restartableRandomGenerator.H"
#include "CarrierField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class turbulentDispersion Declaration
\*---------------------------------------------------------------------------*/

class turbulentDispersion
:
    public LagrangianModel,
    private cloudLagrangianModel
{
private:

    // Private Member Data

        //- Pointer to the drag model
        const drag* dragPtr_;

        //- Reference to the momentum transport model
        const momentumTransportModel& momentumTransportModel_;

        //- Carrier turbulent kinetic energy
        const CarrierField<scalar>& kc_;

        //- Carrier turbulent kinetic energy dissipation rate
        const CarrierField<scalar>& epsilonc_;

        //- Fraction completed of the current turbulent eddy
        Tuple2<bool, CloudStateField<scalar>&> fractionTurb_;

        //- Turbulent time-scale
        Tuple2<bool, CloudStateField<scalar>&> tTurb_;

        //- Turbulent velocity fluctuation
        Tuple2<bool, CloudStateField<vector>&> Uturb_;

        //- Random generator
        restartableRandomGenerator rndGen_;

        //- Average turbulent velocity fluctuation over the current sub-step
        autoPtr<LagrangianSubVectorField> avgUturbPtr_;


    // Private Member Functions

        //- Initialise and return a reference to one of the "turb" fields
        template<class Type>
        Tuple2<bool, CloudStateField<Type>&> initialiseTurbField
        (
            const word& name,
            const dimensionSet& dims,
            const Type& value
        );

        //- Complete the "turb" field, adding any source conditions that are
        //  needed for injections. This has to happen after construction as it
        //  depends on the list of all models.
        template<class InjectionFieldSourceType, class Type>
        void completeTurbField(Tuple2<bool, CloudStateField<Type>&>&);


public:

    //- Runtime type information
    TypeName("turbulentDispersion");


    // Constructors

        //- Construct from components
        turbulentDispersion
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Return the name of the velocity field
        virtual wordList addSupFields() const;

        //- Do post construction steps which require access to other models
        virtual void postConstruct();

        //- Solve equations and/or update continually changing properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Add a source term to the velocity equation
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubVectorSubField& U,
            LagrangianEqn<vector>& eqn
        ) const;

        //- Add a source term to the momentum equation
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& vOrM,
            const LagrangianSubVectorSubField& U,
            LagrangianEqn<vector>& eqn
        ) const;

        //- Write state
        virtual void writeProcessorState(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
