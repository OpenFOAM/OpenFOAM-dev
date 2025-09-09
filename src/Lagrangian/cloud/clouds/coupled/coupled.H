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
    Foam::clouds::coupled

Description
    Base class for clouds which are coupled to a carrier

SourceFiles
    coupled.C

\*---------------------------------------------------------------------------*/

#ifndef coupled_H
#define coupled_H

#include "carried.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class CarrierEqn;

namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class coupled Declaration
\*---------------------------------------------------------------------------*/

class coupled
:
    public carried
{
private:

    // Private Data

        //- Reference to the cloud
        const cloud& cloud_;

        //- Cached carrier kinematic viscosity field
        tmp<volScalarField> tnucVf_;

        //- Carrier equations
        #define DECLARE_CARRIER_EQNS(Type, nullArg)    \
            mutable HashPtrTable<CarrierEqn<Type>>     \
                CAT3(carrier, CAPITALIZE(Type), Eqns_);
        FOR_ALL_FIELD_TYPES(DECLARE_CARRIER_EQNS);
        #undef DECLARE_CARRIER_EQNS

        //- Carrier velocity time derivative
        mutable autoPtr<volVectorField> dUdtcPtr_;


    // Private Member Functions

        //- Get the carrier kinematic viscosity field
        tmp<volScalarField> getNucVf() const;

        //- Calculate the carrier kinematic viscosity
        virtual tmp<LagrangianSubScalarField> calcNuc
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;

        //- Add/get a carrier equation to/from the cache
        template<class Type>
        CarrierEqn<Type>& carrierEqn
        (
            const CloudDerivedField<Type>&,
            const LagrangianSubMesh&
        );

        //- Add/get a carrier equation to/from the cache
        template<class Type>
        CarrierEqn<Type>& carrierEqn
        (
            const CarrierField<Type>&,
            const LagrangianSubMesh&
        );

        //- Access the carrier equations. Only specialisations are defined.
        template<class Type>
        HashPtrTable<CarrierEqn<Type>>& carrierEqns() const;


protected:

    // Protected Member Functions

        //- Update the cached carrier fields
        void updateCarrier();

        //- Clear the carrier equations
        void clearCarrierEqns();

        //- Do re-calculation of a carrier equation for modified particles
        template<class Fieldc>
        bool reCalculateModified
        (
            const LagrangianSubScalarSubField& vOrM,
            const Fieldc& fieldc
        );

        //- Do calculation of a carrier equation
        template<class Fieldc, class Scale = one>
        void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final,
            const LagrangianSubScalarSubField& vOrM,
            const Fieldc& fieldc,
            const Scale& scale = one()
        );

        //- Do calculation of a carrier equation
        template<class Type, class Fieldc, class Scale = one>
        void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final,
            const LagrangianSubScalarSubField& vOrM,
            const LagrangianSubSubField<Type>& field,
            const Fieldc& fieldc,
            const Scale& scale = one()
        );


public:

    // Public Static Data

        //- Run-time type information
        TypeName("coupled");


    // Public Data

        //- Carrier kinematic viscosity
        const CloudDerivedField<scalar>& nuc;


    // Constructors

        //- Construct from a reference to the cloud
        coupled(const cloud& c, const dictionary& dict);


    //- Destructor
    virtual ~coupled();


    // Member Functions

        // Equations

            //- Return if a carrier equation is in the cache
            template<class Type>
            bool hasCarrierEqn(const word&) const;

            //- Get a carrier equation from the cache
            template<class Type>
            const CarrierEqn<Type>& carrierEqn(const word&) const;

            //- Get multiple phase carrier equations from the cache
            template<class Type>
            const HashTable<const CarrierEqn<Type>*> carrierEqns
            (
                const word&
            ) const;

            //- Return if a carrier equation is in the cache
            template<class Type>
            bool hasCarrierEqn(const VolField<Type>&) const;

            //- Get a carrier equation from the cache
            template<class Type>
            const CarrierEqn<Type>& carrierEqn(const VolField<Type>&) const;

            //- Get multiple phase carrier equations from the cache
            template<class Type>
            const HashTable<const CarrierEqn<Type>*> carrierEqns
            (
                const VolField<Type>&
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coupledTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
