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
    Foam::multicomponentCarrierThermo

Description
    Base-class for multicomponent thermodynamic models of the carrier phase of
    a Lagrangian cloud

SourceFiles
    multicomponentCarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentCarrierThermo_H
#define multicomponentCarrierThermo_H

#include "carrierThermo.H"
#include "multicomponentLagrangianSubThermo.H"
#include "speciesTable.H"
#include "DimensionedFieldListSlicer.H"
#include "CarrierThermo.H"
#include "MulticomponentCarrierThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidMulticomponentThermo;

/*---------------------------------------------------------------------------*\
                   Class multicomponentCarrierThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentCarrierThermo
:
    virtual public carrierThermo,
    virtual public multicomponentLagrangianSubThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("multicomponentCarrierThermo");


    //- Destructor
    virtual ~multicomponentCarrierThermo();


    // Member Functions

        // Mass Fractions

            //- Mass fractions []
            virtual const UPtrList<const CarrierField<scalar>>& Y() const = 0;

            //- Mass-fraction for a sub-mesh []
            virtual inline tmp<LagrangianSubScalarSubField> Y
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const;

            //- Mass-fraction for a sub-mesh []
            virtual inline tmp<LagrangianSubScalarSubField> Y
            (
                const word& specieName,
                const LagrangianSubMesh& subMesh
            ) const;
};


/*---------------------------------------------------------------------------*\
             Class multicomponentCarrierThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class multicomponentCarrierThermo::implementation
:
    virtual public multicomponentCarrierThermo
{
    // Private Data

        //- Eulerian thermo
        const fluidMulticomponentThermo& thermo_;

        //- Mass fractions
        UPtrList<const CarrierField<scalar>> Y_;


public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        implementation
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        );

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Species Set

            //- Return the table of species
            virtual const speciesTable& species() const;


        // Mass Fractions

            //- Mass fractions []
            virtual const UPtrList<const CarrierField<scalar>>& Y() const;


        // Composition Access

            //- Get the slicer
            virtual inline DimensionedFieldListSlicer<scalar> Yslicer
            (
                const LagrangianSubMesh& subMesh
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentCarrierThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
