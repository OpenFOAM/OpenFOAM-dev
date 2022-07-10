/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::compressible::thermalBaffle1DFvPatchScalarField

Description
    This BC solves a steady 1D thermal baffle.

    The solid properties are specify as dictionary. Optionally radiative heat
    flux (qr) can be incorporated into the balance. Some under-relaxation might
    be needed on qr.  Baffle and solid properties need to be specified on the
    master side of the baffle.

Usage
    Example of the boundary condition specification using constant
    solid thermo :

    \verbatim
    <masterPatchName>
    {
        type   compressible::thermalBaffle1D<eConstSolidThermoPhysics>;

        samplePatch     <slavePatchName>;

        thickness       uniform 0.005;  // Thickness [m]
        qs              uniform 100;    // Source heat flux [W/m^2]

        qr              none;
        qrRelaxation    1;

        // Solid thermo
        specie
        {
            molWeight       20;
        }
        transport
        {
            kappa           1;
        }
        thermodynamics
        {
            Hf              0;
            Cv              10;
        }
        equationOfState
        {
            rho             10;
        }

        value           uniform 300;
    }

    <slavePatchName>
    {
        type   compressible::thermalBaffle1D<eConstSolidThermoPhysics>;

        samplePatch     <masterPatchName>;

        qr              none;
        relaxation      1;

        value           uniform 300;
    }
    \endverbatim

SourceFiles
    thermalBaffle1DFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffle1DFvPatchScalarField_H
#define thermalBaffle1DFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "autoPtr.H"

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
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of the temperature field
        word TName_;

        //- Baffle is activated
        bool baffleActivated_;

        //- Baffle thickness [m]
        mutable scalarField thickness_;

        //- Source heat flux [W/m^2]
        mutable scalarField qs_;

        //- Solid dictionary
        dictionary solidDict_;

        //- Solid thermo
        mutable autoPtr<solidType> solidPtr_;

        //- Cache qr for relaxation
        scalarField qrPrevious_;

        //- Relaxation for qr
        scalar qrRelaxation_;

        //- Name of the radiative heat flux in local region
        const word qrName_;


    // Private members

        //- Return const solid thermo
        const solidType& solid() const;

        //- Return source heat flux from master
        tmp<scalarField> qs() const;

        //- Return thickness from master
        tmp<scalarField> baffleThickness() const;

        //- Is this patch the owner?
        bool owner() const;

        //- Access the neighbour field
        const thermalBaffle1DFvPatchScalarField& nbrField() const;


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

        //- Disallow copy without setting internal field reference
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
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


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


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
    #include "thermalBaffle1DFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
