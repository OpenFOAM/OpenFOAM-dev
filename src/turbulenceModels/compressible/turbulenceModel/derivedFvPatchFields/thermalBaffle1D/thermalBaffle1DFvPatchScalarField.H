/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::thermalBaffle1DFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description

    This BC solves a steady 1D thermal baffle. The solid properties are
    specify as dictionary. Optionaly radiative heat flux (Qr) can be
    incorporated into the balance. Some under-relaxation might be needed on
    Qr.
    Baffle and solid properties need to be specified on the master side
    of the baffle.

    \heading Patch usage
    Example of the boundary condition specification using constant
    solid thermo :

    \verbatim
    myPatch_master
    {
        type   compressible::thermalBaffle1D<hConstSolidThermoPhysics>;
        samplePatch     myPatch_slave;

        thickness       uniform 0.005;  // thickness [m]
        Qs              uniform 100;    // heat flux [W/m2]

        Qr              none;
        relaxation      0;

        // Solid thermo
        specie
        {
            nMoles          1;
            molWeight       20;
        }
        transport
        {
            kappa           1;
        }
        thermodynamics
        {
            Hf              0;
            Cp              10;
        }
        equationOfState
        {
            rho             10;
        }

        value               uniform 300;
    }

    myPatch_slave
    {
        type   compressible::thermalBaffle1D<hConstSolidThermoPhysics>;
        samplePatch     myPatch_master_master;

        Qr              none;
        relaxation      0;
    }
    \endverbatim



SourceFiles
    thermalBaffle1DFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffle1DFvPatchScalarField_H
#define thermalBaffle1DFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "autoPtr.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class thermalBaffle1DFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

template<class solidType>
class thermalBaffle1DFvPatchScalarField
:
    public mappedPatchBase,
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Baffle is activated
        bool baffleActivated_;

        //- Baffle thickness [m]
        mutable scalarField thickness_;

        //- Superficial heat source [W/m2]
        mutable scalarField Qs_;

        //- Solid dictionary
        dictionary  solidDict_;

        //- Solid thermo
        mutable autoPtr<solidType> solidPtr_;

        //- Chache Qr for relaxation
        scalarField QrPrevious_;

        //- Relaxation for Qr
        scalar QrRelaxation_;

        //- Name of the radiative heat flux in local region
        const word QrName_;


    // Private members

        //- Return const solid thermo
        const solidType& solid() const;

        //- Return Qs from master
        tmp<scalarField> Qs() const;

        //- Return thickness from master
        tmp<scalarField> baffleThickness() const;

        //- Is Owner
        bool owner() const;


public:

    //- Runtime type information
    TypeName("compressible::thermalBaffle1D");


    // Constructors

        //- Construct from patch and internal field
        thermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermalBaffle1DFvPatchScalarField onto a new patch
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalBaffle1DFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&,
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
                new thermalBaffle1DFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "thermalBaffle1DFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
