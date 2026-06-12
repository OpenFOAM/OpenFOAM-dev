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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;
class multicomponentThermo;
class fluidThermo;
class fluidMulticomponentThermo;

namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class coupledToThermalFluid Declaration
\*---------------------------------------------------------------------------*/

class coupledToThermalFluid
:
    public coupledToFluid
{
private:

    // Private Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;

        //- Carrier thermo
        const fluidThermo& thermoc_;

        //- Carrier multicomponent thermo (may be null)
        const fluidMulticomponentThermo& multicomponentThermoc_;

        //- Corresponding Eulerian thermo (may be null)
        const basicThermo& thermocPhase_;

        //- Corresponding multicomponent Eulerian thermo (may be null)
        const multicomponentThermo& multicomponentThermocPhase_;


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

        //- Corresponding Eulerian enthalpy/energy
        const CarrierField<scalar>& hecPhase;

        //- Carrier heat capacity
        const CarrierField<scalar>& Cpvc;

        //- Carrier heat capacity at constant pressure
        const CarrierField<scalar>& Cpc;

        //- Carrier heat capacity at constant volume
        const CarrierField<scalar>& Cvc;

        //- Carrier thermal conductivity
        const CarrierField<scalar>& kappac;

        //- Carrier Prandtl number
        const CloudDerivedField<scalar>& Prc;

        //- Map from cloud to carrier specie indices
        labelList iToic;

        //- Map from cloud to corresponding Eulerian specie indices
        labelList iToicPhase;

        //- Carrier mass fractions (indexed with the carrier specie index)
        UPtrList<const CarrierField<scalar>> Yc;

        //- Corresponding Eulerian mass fractions (indexed with the Eulerian
        //  phase specie index)
        UPtrList<const CarrierField<scalar>> YcPhase;


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        coupledToThermalFluid(const cloud& c, const carried&, const thermal&);


    //- Destructor
    virtual ~coupledToThermalFluid();


    // Member Functions

        //- Return whether this cloud has a carrier thermo
        bool hasThermoc() const;

        //- Return whether this cloud has a corresponding Eulerian phase thermo
        bool hasThermocPhase() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
