/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "tractionDisplacementFvPatchVectorField.H"
#include "solidDisplacementThermo.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::tractionDisplacementFvPatchVectorField::updateCoeffs
(
    const Type& pressure
)
{
    const label patchi = patch().index();

    const solidDisplacementThermo& thermo =
        db().lookupObject<solidDisplacementThermo>
        (
            physicalProperties::typeName
        );

    const scalarField& E = thermo.E(patchi);
    const scalarField& nu = thermo.nu(patchi);

    const scalarField mu(E/(2.0*(1.0 + nu)));
    const scalarField lambda
    (
        thermo.planeStress()
      ? nu*E/((1 + nu)*(1 - nu))
      : nu*E/((1 + nu)*(1 - 2*nu))
    );
    const scalarField threeK
    (
        thermo.planeStress()
      ? E/(1 - nu)
      : E/(1 - 2*nu)
    );

    const scalarField twoMuLambda(2*mu + lambda);

    const vectorField n(patch().nf());

    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    gradient() =
    (
        (traction_ - pressure*n)
      + twoMuLambda*fvPatchField<vector>::snGrad() - (n & sigmaD)
    )/twoMuLambda;

    if (thermo.thermalStress())
    {
        const scalarField& alphav = thermo.alphav(patchi);

        gradient() +=
            n*threeK*alphav*thermo.T().boundaryField()[patchi]/twoMuLambda;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// ************************************************************************* //
