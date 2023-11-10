/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::multiphaseExternalTemperatureFvPatchScalarField

Description
    Mixed boundary condition for the phase temperature of a phase in an
    Euler-Euler multiphase simulation, to be used for heat-transfer with another
    region in a CHT case.  Optional thin wall material layer resistances can be
    specified through thicknessLayers and kappaLayers entries.

See also
    Foam::externalTemperatureFvPatchScalarField

SourceFiles
    multiphaseExternalTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseExternalTemperatureFvPatchScalarField_H
#define multiphaseExternalTemperatureFvPatchScalarField_H

#include "externalTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class multiphaseExternalTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class multiphaseExternalTemperatureFvPatchScalarField
:
    public externalTemperatureFvPatchScalarField
{

protected:

    //- Get the patch kappa, sum kappa*Tc/delta, kappa/delta and
    //  reference temperature for all phases except the phase being solved
    //  and also the current phase-average wall T and phase heat-flux
    //  obtained by partitioning the sum heat-flux provided
    virtual void getKappa
    (
        scalarField& kappa,
        scalarField& sumKappaTByDelta,
        scalarField& sumKappaByDelta,
        scalarField& Tref,
        scalarField& Tw,
        scalarField& sumq,
        scalarField& qByKappa
    ) const;


public:

    //- Runtime type information
    TypeName("multiphaseExternalTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        multiphaseExternalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  multiphaseExternalTemperatureFvPatchScalarField
        //  onto a new patch
        multiphaseExternalTemperatureFvPatchScalarField
        (
            const multiphaseExternalTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        multiphaseExternalTemperatureFvPatchScalarField
        (
            const multiphaseExternalTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        multiphaseExternalTemperatureFvPatchScalarField
        (
            const multiphaseExternalTemperatureFvPatchScalarField&,
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
                new multiphaseExternalTemperatureFvPatchScalarField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
