/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::multiphaseCoupledTemperatureFvPatchScalarField

Description
    Mixed boundary condition for the phase temperature of a phase in an
    Euler-Euler multiphase simulation, to be used for heat-transfer with another
    region in a CHT case.  Optional thin wall material layer resistances can be
    specified through thicknessLayers and kappaLayers entries.

See also
    Foam::coupledTemperatureFvPatchScalarField

SourceFiles
    multiphaseCoupledTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseCoupledTemperatureFvPatchScalarField_H
#define multiphaseCoupledTemperatureFvPatchScalarField_H

#include "coupledTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class multiphaseCoupledTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class multiphaseCoupledTemperatureFvPatchScalarField
:
    public coupledTemperatureFvPatchScalarField
{

protected:

    // Protected Member Functions

        //- Get the patch kappa, sum kappa*Tc/delta and kappa/delta for all
        //  phases except the phase being solved and also the phase heat-flux
        //  obtained by partitioning the sum heat-flux provided
        virtual void getThis
        (
            tmp<scalarField>& kappa,
            tmp<scalarField>& sumKappaTByDelta,
            tmp<scalarField>& sumKappaByDelta,
            scalarField& sumq,
            tmp<scalarField>& qByKappa
        ) const;

        //- Get the neighbour patch sum kappa*Tc/delta and kappa/delta
        //  for all phases and the heat-flux correction
        virtual void getNbr
        (
            tmp<scalarField>& sumKappaTByDeltaNbr,
            tmp<scalarField>& sumKappaByDeltaNbr,
            tmp<scalarField>& qNbr
        ) const;

        //- Get the neighbour patch wall temperature and heat-flux
        //  when wall thermal resistance is specified
        virtual void getNbr
        (
            tmp<scalarField>& TwNbr,
            tmp<scalarField>& qNbr
        ) const;


public:

    //- Runtime type information
    TypeName("multiphaseCoupledTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        multiphaseCoupledTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  multiphaseCoupledTemperatureFvPatchScalarField
        //  onto a new patch
        multiphaseCoupledTemperatureFvPatchScalarField
        (
            const multiphaseCoupledTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        multiphaseCoupledTemperatureFvPatchScalarField
        (
            const multiphaseCoupledTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        multiphaseCoupledTemperatureFvPatchScalarField
        (
            const multiphaseCoupledTemperatureFvPatchScalarField&,
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
                new multiphaseCoupledTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Operators

        //- Inherit assignment
        using coupledTemperatureFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
