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
    Foam::multicomponentSurfaceThermo

Description
    Base-class for multicomponent Lagrangian surface thermodynamic models

SourceFiles
    multicomponentSurfaceThermoI.H
    multicomponentSurfaceThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentSurfaceThermo_H
#define multicomponentSurfaceThermo_H

#include "surfaceThermo.H"
#include "multicomponentLagrangianSubThermo.H"
#include "speciesTable.H"
#include "DimensionedFieldListSlicer.H"
#include "SurfaceThermo.H"
#include "MulticomponentSurfaceThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidMulticomponentThermo;

/*---------------------------------------------------------------------------*\
                   Class multicomponentSurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentSurfaceThermo
:
    virtual public surfaceThermo,
    virtual public multicomponentLagrangianSubThermo
{
protected:

    // Protected Member Functions

        //- Reset to the far thermo
        virtual void reset() const;

        //- Reset to the far thermo
        virtual void resetYs() const = 0;


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("multicomponentSurfaceThermo");


    //- Destructor
    virtual ~multicomponentSurfaceThermo();


    // Member Functions

        // Species Set

            //- Return the table of species
            virtual const speciesTable& species() const = 0;

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const;


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

            //- Allow a model to set the mass fractions []
            virtual void setYs
            (
                const LagrangianModel& model,
                PtrList<LagrangianSubScalarSubField>&& Ys
            ) = 0;

            //- Allow a model to set the mass fractions []
            virtual void setYs
            (
                const LagrangianModel& model,
                PtrList<LagrangianSubScalarField>&& Ys
            );
};


/*---------------------------------------------------------------------------*\
             Class multicomponentSurfaceThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class multicomponentSurfaceThermo::implementation
:
    virtual public multicomponentSurfaceThermo
{
protected:

    // Protected Data

        //- Bulk thermo
        const multicomponentLagrangianSubThermo& farMulticomponentThermo_;

        //- Surface mass fractions
        mutable PtrList<LagrangianSubScalarSubField> Ys_;

        //- The model that set the surface mass fractions
        mutable const LagrangianModel* YsModelPtr_;


    // Protected Member Functions

        //- Reset to the far thermo
        virtual void resetYs() const;


public:

    // Constructors

        //- Construct from a reference to the far thermo
        implementation(const basicLagrangianSubThermo& farThermo);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Species Set

            //- Return the table of species
            virtual const speciesTable& species() const;


        // Mass Fractions

            //- Mass-fraction for a sub-mesh []
            virtual tmp<LagrangianSubScalarSubField> Y
            (
                const label speciei,
                const LagrangianSubMesh& subMesh
            ) const;

            //- Mass-fraction for a sub-mesh []
            virtual tmp<LagrangianSubScalarSubField> Y
            (
                const word& specieName,
                const LagrangianSubMesh& subMesh
            ) const;

            //- Allow a model to set the mass fractions []
            virtual void setYs
            (
                const LagrangianModel& model,
                PtrList<LagrangianSubScalarSubField>&& Ys
            );


        // Composition Access

            //- Get the slicer
            inline DimensionedFieldListSlicer<scalar> Yslicer
            (
                const LagrangianSubMesh& subMesh
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentSurfaceThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
