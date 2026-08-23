/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::clouds::coupledToThermalFluid

Description
    Base class for clouds which are coupled to a fluid with a thermodynamic
    model

SourceFiles
    coupledToThermalFluid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledToThermalFluid_H
#define coupledToThermalFluid_H

#include "thermal.H"
#include "coupledToFluid.H"
#include "fluidCarrierThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                     Class coupledToThermalFluid Declaration
\*---------------------------------------------------------------------------*/

class coupledToThermalFluid
:
    private autoPtr<fluidCarrierThermo>,
    private autoPtr<fluidSurfaceThermo>,
    private autoPtr<carrierThermo>,
    public coupledToFluid
{
private:

    // Private Static Data

        //- Null/empty set of mass fractions
        static const UPtrList<const CarrierField<scalar>> nullY_;


    // Private Data

        //- Reference to the cloud
        const cloud& cloud_;


protected:

    // Protected Member Functions

        //- Update the cached carrier fields
        void updateCarrier();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("coupledToThermalFluid");


    // Public Data

        //- Carrier pressure
        const CarrierField<scalar>& pc;

        //- Carrier temperature
        const CarrierField<scalar>& Tc;

        //- Carrier enthalpy/energy
        const CarrierField<scalar>& hec;

        //- Corresponding Eulerian phase enthalpy/energy
        const CarrierField<scalar>& hecPhase;

        //- Carrier mass fractions (indexed with the carrier specie index)
        const UPtrList<const CarrierField<scalar>>& Yc;

        //- Corresponding Eulerian phase mass fractions
        //  (indexed with the Eulerian phase specie index)
        const UPtrList<const CarrierField<scalar>>& YcPhase;

        //- Map from cloud to carrier specie indices
        labelList iToic;

        //- Map from cloud to corresponding Eulerian phase specie indices
        labelList iToicPhase;


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        coupledToThermalFluid(const cloud& c, const carried&, const thermal&);


    //- Destructor
    virtual ~coupledToThermalFluid();


    // Member Functions

        //- Return whether this cloud has a carrier thermo
        bool hasThermoc() const;

        //- Return whether this cloud has a carrier surface thermo
        bool hasThermocs() const;

        //- Return whether this cloud has a corresponding Eulerian phase thermo
        bool hasThermocPhase() const;

        //- Determine whether the carrier thermo model is one of the given types
        template<class Thermo, class ... Thermos>
        bool isThermoc() const;

        //- Determine whether the carrier surface thermo model is one of the
        //  given types
        template<class Thermo, class ... Thermos>
        bool isThermocs() const;

        //- Determine whether the corresponding Eulerian phase thermo model is
        //  one of the given types
        template<class Thermo, class ... Thermos>
        bool isThermocPhase() const;

        //- Generate an error if the carrier thermo model is not one of the
        //  given types
        template<class Thermo, class ... Thermos>
        void assertThermoc() const;

        //- Generate an error if the carrier surface thermo model is not one of
        //  the given types
        template<class Thermo, class ... Thermos>
        void assertThermocs() const;

        //- Generate an error if the corresponding Eulerian phase thermo model
        //  is not one of the given types
        template<class Thermo, class ... Thermos>
        void assertThermocPhase() const;

        //- Generate an error if the carrier thermo model is not one of the
        //  given types. A model is provided to give the error message some
        //  context.
        template<class Thermo, class ... Thermos>
        void assertThermoc(const LagrangianModel& model) const;

        //- Generate an error if the carrier surface thermo model is not one of
        //  the given types. A model is provided to give the error message some
        //  context.
        template<class Thermo, class ... Thermos>
        void assertThermocs(const LagrangianModel& model) const;

        //- Generate an error if the corresponding Eulerian phase thermo model
        //  is not one of the given types. A model is provided to give the
        //  error message some context.
        template<class Thermo, class ... Thermos>
        void assertThermocPhase(const LagrangianModel& model) const;

        //- Get a reference to the carrier thermo model of a given type
        template<class Thermo = fluidCarrierThermo, class ... Args>
        const Thermo& thermoc(const Args& ... args) const;

        //- Get a reference to the carrier surface thermo model of a given type
        template<class Thermo = fluidSurfaceThermo, class ... Args>
        const Thermo& thermocs(const Args& ... args) const;

        //- Get a non-const reference to the carrier surface thermo model of a
        //  given type
        template<class Thermo = fluidSurfaceThermo, class ... Args>
        Thermo& thermocsRef(const Args& ... args);

        //- Get a reference to the corresponding Eulerian phase thermo model of
        //  a given type
        template<class Thermo = carrierThermo, class ... Args>
        const Thermo& thermocPhase(const Args& ... args) const;

        //- Return the carrier thermal conductivity
        tmp<LagrangianSubScalarSubField> kappac(const LagrangianSubMesh&) const;

        //- Return the carrier Prandtl number
        tmp<LagrangianSubScalarField> Prc(const LagrangianSubMesh&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coupledToThermalFluidTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
