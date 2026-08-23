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
    Foam::basicLagrangianSubThermo

Description
    Abstract interface which provides access to basic (i.e., applicable to both
    solids and fluids) thermodynamic properties on Lagrangian sub-meshes.

SourceFiles
    basicLagrangianSubThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicLagrangianSubThermo_H
#define basicLagrangianSubThermo_H

#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOdictionary;
class LagrangianMesh;

/*---------------------------------------------------------------------------*\
                   Class basicLagrangianSubThermo Declaration
\*---------------------------------------------------------------------------*/

class basicLagrangianSubThermo
{
public:

    //- Runtime type information
    TypeName("basicLagrangianSubThermo");


    //- Destructor
    virtual ~basicLagrangianSubThermo();


    // Member Functions

        //- Access the properties dictionary
        virtual const IOdictionary& properties() const = 0;

        //- Access the mesh
        virtual const LagrangianMesh& mesh() const = 0;

        //- Phase name
        virtual const word& phaseName() const = 0;

        //- Return the name of the mixture
        virtual word mixtureName() const = 0;


        // Molecular Properties

            //- Molecular weight for a sub-mesh [kg/kmol]
            virtual tmp<LagrangianSubScalarField> W
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;


        // Thermodynamic State

            //- Temperature for a sub-mesh [K]
            virtual tmp<LagrangianSubScalarSubField> T
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;


        // Thermodynamic Properties

            //- Density for a sub-mesh [kg/m^3]
            virtual tmp<LagrangianSubScalarSubField> rho
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Internal energy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarSubField> e
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Sensible enthalpy for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> hs
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Absolute enthalpy for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> ha
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Heat capacity at constant pressure for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cp
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Heat capacity at constant volume for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarSubField> Cv
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Coefficient of thermal expansion for a sub-mesh [1/K]
            virtual tmp<LagrangianSubScalarField> alphav
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;


        // Transport Properties

            //- Thermal conductivity for a sub-mesh [W/m/K]
            virtual tmp<LagrangianSubScalarSubField> kappa
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Thermal dynamic diffusivity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarField> kappaByCp
            (
                const LagrangianSubMesh& subMesh
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
