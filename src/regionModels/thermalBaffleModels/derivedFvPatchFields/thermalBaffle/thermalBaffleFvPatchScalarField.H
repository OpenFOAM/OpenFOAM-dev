/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::compressible::thermalBaffleFvPatchScalarField

Description
    This boundary condition provides a coupled temperature condition between
    multiple mesh regions.

    The regions are generally referred to as the:
      - primary region,
      - and baffle region.

    The primary region creates the baffle region and evolves its energy
    equation either:
      - 1-D, normal to each patch face
      - 2-D, normal and tangential components

    The thermodynamic properties of the baffle material are specified via
    dictionary entries on the master patch.

Usage
    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        type                compressible::thermalBaffle;

        // Underlaying coupled boundary condition
        Tnbr               T;
        kappa              fluidThermo; // or solidThermo
        KappaName          none;
        qrNbr              qr;//or none.Name of qr field on neighbourregion
        qr                 none;// or none.Name of qr field on localregion
        value              uniform 300;

        // Baffle region name
        regionName          baffleRegion;
        active              yes;

        // Solid thermo in solid region
        thermoType
        {
            type            heSolidThermo;
            mixture         pureMixture;
            transport       constIso;
            thermo          hConst;
            equationOfState rhoConst;
            specie          specie;
            energy          sensibleEnthalpy;
        }

        mixture
        {
            specie
            {
                molWeight       20;
            }
            transport
            {
                kappa           0.01;
            }
            thermodynamics
            {
                Hf              0;
                Cp              15;
            }
            density
            {
                rho             80;
            }
        }

        radiation
        {
            radiationModel  opaqueSolid;
            absorptionEmissionModel none;
            scatterModel    none;
        }

        // Extrude model for new region

        extrudeModel        linearNormal;
        nLayers             50;
        expansionRatio      1;
        columnCells         false; //3D or 1D
        linearNormalCoeffs
        {
            thickness           0.02;
        }

    }
    \endverbatim

    Slave patch on primary region:
    \verbatim
    <slavePatchName>
    {
        type                compressible::thermalBaffle;
        kappa               fluidThermo;
        kappaName           none;
        value               uniform 300;
    \endverbatim

    Patches on baffle region:
    \verbatim
    bottom
    {
        type                compressible::thermalBaffle;
        kappa               solidThermo;
        kappaName           none;
        value               uniform 300;
    }

    top
    {
        type                compressible::thermalBaffle;
        kappa               solidThermo;
        kappaName           none;
        value               uniform 300;
    }
    \endverbatim

See also
    Foam::compressible::turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
    Foam::regionModels::thermalBaffleModels::thermalBaffleModel

SourceFiles
    thermalBaffleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffleFvPatchScalarField_H
#define thermalBaffleFvPatchScalarField_H


#include "autoPtr.H"
#include "regionModel.H"
#include "thermalBaffleModel.H"
#include "extrudePatchMesh.H"
#include "turbulentTemperatureRadCoupledMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
        Class thermalBaffleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalBaffleFvPatchScalarField
:
    public turbulentTemperatureRadCoupledMixedFvPatchScalarField
{
    // Private Data

        //- Enumeration of patch IDs
        enum patchID
        {
            bottomPatchID,
            topPatchID,
            sidePatchID
        };

        //- Is the baffle owner
        bool owner_;

        //- Thermal baffle
        autoPtr<regionModels::thermalBaffleModels::thermalBaffleModel> baffle_;

        //- Dictionary
        dictionary dict_;

        //- Auto pointer to extrapolated mesh from patch
        autoPtr<extrudePatchMesh> extrudeMeshPtr_;


    // Private Member Functions

        //- Extrude mesh
        void createPatchMesh();

public:

    //- Runtime type information
    TypeName("compressible::thermalBaffle");


    // Constructors

        //- Construct from patch and internal field
        thermalBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermalBaffleFvPatchScalarField onto a new patch
        thermalBaffleFvPatchScalarField
        (
            const thermalBaffleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalBaffleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalBaffleFvPatchScalarField
        (
            const thermalBaffleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalBaffleFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
