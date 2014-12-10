/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::energyRegionCoupledFvPatchScalarField

Description
    Energy region coupled implicit boundary condition.
    The fvPatch is treated as uncoupled from the delta point of view.
    In the mesh the fvPatch is an interface and is incorporated
    into the matrix implicitly.

SourceFiles
    energyRegionCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyRegionCoupledFvPatchScalarField_H
#define energyRegionCoupledFvPatchScalarField_H

#include "regionCoupledBaseFvPatch.H"
#include "LduInterfaceField.H"
#include "fvPatchField.H"
#include "NamedEnum.H"
#include "basicThermo.H"
#include "coupledFvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class energyRegionCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class energyRegionCoupledFvPatchScalarField
:
    public coupledFvPatchField<scalar>
{

public:

        enum kappaMethodType
        {
            SOLID,
            FLUID,
            UNDEFINED
        };


private:

    // Private data

        //- Local reference to region couple patch
        const regionCoupledBaseFvPatch& regionCoupledPatch_;

        //- Methof to extract kappa
        static const NamedEnum<kappaMethodType, 3> methodTypeNames_;

        //- How to get K
        mutable kappaMethodType method_;

        //- AutoPtr to nbr thermo
        mutable const basicThermo* nbrThermoPtr_;

        //- AutoPtr to my thermo
        mutable const basicThermo* thermoPtr_;


    // Private functions

        //- Local weight for this coupled field
        tmp<scalarField> weights() const;

        //- Return nbr temperature internal field
        tmp<scalarField> patchNeighbourTemperatureField() const;

        //- Return local temperature internal field
        tmp<scalarField> patchInternalTemperatureField() const;

        //- Return kappa
        tmp<scalarField> kappa() const;

        //- Set method
        void setMethod() const;


public:

    //- Runtime type information
    TypeName("compressible::energyRegionCoupled");


    // Constructors

        //- Construct from patch and internal field
        energyRegionCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energyRegionCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given energyRegionCoupledFvPatchScalarField
        // onto a new patch
        energyRegionCoupledFvPatchScalarField
        (
            const energyRegionCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        energyRegionCoupledFvPatchScalarField
        (
            const energyRegionCoupledFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new energyRegionCoupledFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        energyRegionCoupledFvPatchScalarField
        (
            const energyRegionCoupledFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new energyRegionCoupledFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~energyRegionCoupledFvPatchScalarField()
    {}


    // Member functions

        // Access

            //- Method to obtain K
            word kappaMethod() const
            {
                return methodTypeNames_[method_];
            }


        // Evaluation functions

            //- Return neighbour coupled internal cell data
            virtual tmp<scalarField> patchNeighbourField() const;

            //- Return patch-normal gradient
            virtual tmp<scalarField> snGrad() const;

            //- Return patch-normal gradient
            //  Note: the deltaCoeffs supplied are not used
            virtual tmp<scalarField> snGrad
            (
                const scalarField& deltaCoeffs
            ) const;


            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );


        // Coupled interface functionality

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<scalar>& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<scalar>&,
                const Field<scalar>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return regionCoupledPatch_.regionCoupleType();
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
