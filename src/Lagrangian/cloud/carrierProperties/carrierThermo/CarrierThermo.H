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
    Foam::CarrierThermo

Description
    Thermodynamic model of the carrier phase of a Lagrangian cloud.

SourceFiles
    CarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef CarrierThermo_H
#define CarrierThermo_H

#include "carried.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CarrierThermoName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(CarrierThermo);


/*---------------------------------------------------------------------------*\
                        Class CarrierThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType, class BasicThermoType>
class CarrierThermo
:
    public CarrierThermoName,
    public MixtureType,
    public BasicThermoType
{
protected:

    // Protected Member Functions

        //- Return a LagrangianSubScalarField of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<LagrangianSubScalarField> LagrangianSubScalarFieldProperty
        (
            const LagrangianSubMesh& subMesh,
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;


public:

    // Typedefs

        //- Mixture type
        typedef MixtureType mixtureType;

        //- Basic thermo
        typedef BasicThermoType basicThermoType;


    //- Disambiguate debug switch used by derivations
    using CarrierThermoName::debug;


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        CarrierThermo
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        );

        //- Disallow default bitwise copy construction
        CarrierThermo(const CarrierThermo&) = delete;


    //- Destructor
    virtual ~CarrierThermo();


    // Member Functions

        //- Return the name of the mixture
        virtual word mixtureName() const;


        // Molecular Properties

            //- Molecular weight for a sub-mesh [kg/kmol]
            virtual tmp<LagrangianSubScalarField> W
            (
                const LagrangianSubMesh&
            ) const;


        // Thermodynamic Properties

            //- Internal energy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarSubField> e
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Sensible enthalpy for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> hs
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Absolute enthalpy for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> ha
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Coefficient of thermal expansion for a sub-mesh [1/K]
            virtual tmp<LagrangianSubScalarField> alphav
            (
                const LagrangianSubMesh& subMesh
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CarrierThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
