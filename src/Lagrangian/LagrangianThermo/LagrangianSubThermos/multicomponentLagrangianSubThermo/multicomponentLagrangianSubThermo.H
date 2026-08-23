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
    Foam::multicomponentLagrangianSubThermo

Description
    Abstract interface which provides access to basic thermodynamic properties
    of individual species within a multicomponent mixture on Lagrangian
    sub-meshes.

SourceFiles
    multicomponentLagrangianSubThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentLagrangianSubThermo_H
#define multicomponentLagrangianSubThermo_H

#include "basicLagrangianSubThermo.H"
#include "speciesTable.H"
#include "DimensionedFieldListSlicer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class multicomponentLagrangianSubThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentLagrangianSubThermo
:
    virtual public basicLagrangianSubThermo
{
public:

    //- Runtime type information
    TypeName("multicomponentLagrangianSubThermo");


    //- Destructor
    virtual ~multicomponentLagrangianSubThermo();


    // Member Functions

        // Species Set

            //- Return the table of species
            virtual const speciesTable& species() const = 0;

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const
            {
                return species().found(specieName);
            }


        // Specie Molecular Properties

            //- Molecular weight [kg/kmol]
            virtual dimensionedScalar Wi(const label speciei) const = 0;


        // Mass Fractions

            //- Mass-fraction for a sub-mesh []
            virtual tmp<LagrangianSubScalarSubField> Y
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Mass-fraction for a sub-mesh []
            virtual tmp<LagrangianSubScalarSubField> Y
            (
                const word& specieName,
                const LagrangianSubMesh& subMesh
            ) const = 0;


        // Specie Thermodynamic Properties

            //- Density for a sub-mesh [kg/m^3]
            virtual tmp<LagrangianSubScalarField> rhoi
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Sensible enthalpy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarField> hsi
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Absolute enthalpy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarField> hai
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual dimensionedScalar hfi(const label speciei) const = 0;


        // Composition Access

            //- Get the slicer for a sub-mesh
            virtual DimensionedFieldListSlicer<scalar> Yslicer
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Get the composition of an element
            inline scalarFieldListSlice elementComposition
            (
                const DimensionedFieldListSlicer<scalar>& Yslicer,
                const label i
            ) const
            {
                return Yslicer.slice(i);
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
