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
    Foam::MulticompoonentCarrierThermo

Description
    Thermodynamic model of the carrier phase of a Lagrangian cloud.

SourceFiles
    MulticompoonentCarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef MulticomponentCarrierThermo_H
#define MulticomponentCarrierThermo_H

#include "carried.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class MulticomponentCarrierThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class MulticomponentCarrierThermo
:
    public BaseThermo
{
protected:

    // Protected Member Functions

        //- Return a LagrangianSubScalarField of the given property
        template<class Method, class ... Args>
        tmp<LagrangianSubScalarField> LagrangianSubScalarFieldPropertyi
        (
            const LagrangianSubMesh& subMesh,
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const label speciei,
            const Args& ... args
        ) const;


public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        MulticomponentCarrierThermo
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        );


    //- Destructor
    virtual ~MulticomponentCarrierThermo();


    // Member Functions

        // Molecular Properties

            //- Molecular weight [kg/kmol]
            virtual dimensionedScalar Wi(const label speciei) const;


        // Thermodynamic Properties

            //- Density for a sub-mesh [kg/m^3]
            virtual tmp<LagrangianSubScalarField> rhoi
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const;

            //- Sensible enthalpy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarField> hsi
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const;

            //- Absolute enthalpy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarField> hai
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const;

            //- Enthalpy of formation [J/kg]
            virtual dimensionedScalar hfi(const label speciei) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MulticomponentCarrierThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
