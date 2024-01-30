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

\*---------------------------------------------------------------------------*/

#include "filmSurfaceVelocityFvPatchVectorField.H"
#include "compressibleMomentumTransportModel.H"
#include "mappedFvPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmSurfaceVelocityFvPatchVectorField::
filmSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF, dict, false),
    Cs_(dict.lookupOrDefault<scalar>("Cs", 0))
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // If the value entry is not present initialise to zero-gradient
        fvPatchVectorField::operator=(patchInternalField());
    }
}


Foam::filmSurfaceVelocityFvPatchVectorField::
filmSurfaceVelocityFvPatchVectorField
(
    const filmSurfaceVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    Cs_(ptf.Cs_)
{}


Foam::filmSurfaceVelocityFvPatchVectorField::
filmSurfaceVelocityFvPatchVectorField
(
    const filmSurfaceVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    Cs_(ptf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmSurfaceVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get the coupling information from the mappedPatchBase
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& patchNbr = mapper.nbrFvPatch();

    // Neighbour patch face-cell velocity
    const vectorField UpNbr
    (
        patchNbr.lookupPatchField<volVectorField, scalar>("U")
       .patchInternalField()
    );

    // Set the reference value to the neighbouring fluid internal velocity
    refValue() = mapper.fromNeighbour(UpNbr);

    // Remove the normal component of the surface vel
    const vectorField n(patch().nf());
    refValue() -= n*(n & refValue());

    // Lookup the momentum transport model
    const compressibleMomentumTransportModel& transportModel =
        db().lookupType<compressibleMomentumTransportModel>();

    // Lookup the neighbour momentum transport model
    const compressibleMomentumTransportModel& transportModelNbr =
        mapper.nbrMesh().lookupType<compressibleMomentumTransportModel>();

    // Patch laminar dynamic viscosity divided by delta
    const tmp<scalarField> muEffByDelta
    (
        transportModel.rho().boundaryField()[patch().index()]
       *transportModel.nuEff(patch().index())
       *patch().deltaCoeffs()
    );

    if (Cs_ > 0)
    {
        // Get the neighbour patch density
        const tmp<scalarField> rhopNbr
        (
            mapper.fromNeighbour
            (
                transportModelNbr.rho().boundaryField()[patchNbr.index()]
            )
        );

        // Calculate the drag coefficient from the drag constant
        // and the magnitude of the velocity difference
        const scalarField Ds(Cs_*rhopNbr*mag(refValue() - *this));

        // Calculate the value-fraction from the balance between the
        // external fluid drag and internal film stress
        valueFraction() = Ds/(muEffByDelta + Ds);
    }
    else
    {
        // Get the neighbour patch laminar dynamic viscosity divided by delta
        const tmp<scalarField> muEffByDeltaNbr
        (
            mapper.fromNeighbour
            (
                transportModelNbr.rho().boundaryField()[patchNbr.index()]
               *transportModelNbr.nuEff(patchNbr.index())
               *patchNbr.deltaCoeffs()
            )
        );

        // Calculate the value-fraction from the balance between the
        // external fluid and internal film stresses
        valueFraction() = muEffByDeltaNbr()/(muEffByDelta + muEffByDeltaNbr());
    }

    mixedFvPatchField<vector>::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::filmSurfaceVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    if (Cs_ > 0)
    {
        writeEntry(os, "Cs", Cs_);
    }

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       filmSurfaceVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
