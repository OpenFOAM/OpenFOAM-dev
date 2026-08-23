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

\*---------------------------------------------------------------------------*/

#include "surfaceFunction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianSubField<Foam::scalar, PrimitiveField>>
Foam::surfaceFunction::transportProperty
(
    const LagrangianSubMesh& subMesh,
    void (surfaceFunction::*coeffs)
    (
        const LagrangianSubMesh& subMesh,
        autoPtr<dimensionedScalar>&,
        tmp<LagrangianSubScalarField>&
    ) const,
    const Thermo& farThermo,
    const Thermo& surfaceThermo,
    tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
    (
        const LagrangianSubMesh& subMesh
    ) const
) const
{
    autoPtr<dimensionedScalar> fValuePtr;
    tmp<LagrangianSubScalarField> tfField;
    (this->*coeffs)(subMesh, fValuePtr, tfField);

    if (fValuePtr.valid())
    {
        if (fValuePtr().value() == 0) return (farThermo.*property)(subMesh);

        if (fValuePtr().value() == 1) return (surfaceThermo.*property)(subMesh);

        return
            toPrimitiveField<PrimitiveField>
            (
               (1 - fValuePtr().value())*(farThermo.*property)(subMesh)
             + fValuePtr().value()*(surfaceThermo.*property)(subMesh)
            );
    }

    return
        toPrimitiveField<PrimitiveField>
        (
            (scalar(1) - tfField())*(farThermo.*property)(subMesh)
          + tfField()*(surfaceThermo.*property)(subMesh)
        );
}


template<class FarFunction, class SurfaceFunction>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceFunction::transportProperty
(
    const LagrangianSubMesh& subMesh,
    void (surfaceFunction::*coeffs)
    (
        const LagrangianSubMesh& subMesh,
        autoPtr<dimensionedScalar>&,
        tmp<LagrangianSubScalarField>&
    ) const,
    const FarFunction& farFunction,
    const SurfaceFunction& surfaceFunction
) const
{
    autoPtr<dimensionedScalar> fValuePtr;
    tmp<LagrangianSubScalarField> tfField;
    (this->*coeffs)(subMesh, fValuePtr, tfField);

    if (fValuePtr.valid())
    {
        if (fValuePtr().value() == 0) return farFunction();

        if (fValuePtr().value() == 1) return surfaceFunction();

        return
            (1 - fValuePtr().value())*farFunction()
          + fValuePtr().value()*surfaceFunction();
    }

    return
        (scalar(1) - tfField())*farFunction()
      + tfField()*surfaceFunction();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianSubField<Foam::scalar, PrimitiveField>>
Foam::surfaceFunction::thermalConductivity
(
    const LagrangianSubMesh& subMesh,
    const Thermo& farThermo,
    const Thermo& surfaceThermo,
    tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
    (
        const LagrangianSubMesh& subMesh
    ) const
) const
{
    return
        transportProperty
        (
            subMesh,
            &surfaceFunction::thermalConductivityCoeffs,
            farThermo,
            surfaceThermo,
            property
        );
}


template<class Thermo, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianSubField<Foam::scalar, PrimitiveField>>
Foam::surfaceFunction::dynamicDiffusivity
(
    const LagrangianSubMesh& subMesh,
    const Thermo& farThermo,
    const Thermo& surfaceThermo,
    tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
    (
        const LagrangianSubMesh& subMesh
    ) const
) const
{
    return
        transportProperty
        (
            subMesh,
            &surfaceFunction::dynamicDiffusivityCoeffs,
            farThermo,
            surfaceThermo,
            property
        );
}


template<class Thermo, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianSubField<Foam::scalar, PrimitiveField>>
Foam::surfaceFunction::kinematicDiffusivity
(
    const LagrangianSubMesh& subMesh,
    const Thermo& farThermo,
    const Thermo& surfaceThermo,
    tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
    (
        const LagrangianSubMesh& subMesh
    ) const
) const
{
    return
        transportProperty
        (
            subMesh,
            &surfaceFunction::kinematicDiffusivityCoeffs,
            farThermo,
            surfaceThermo,
            property
        );
}


template<class FarFunction, class SurfaceFunction>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceFunction::kinematicDiffusivity
(
    const LagrangianSubMesh& subMesh,
    const FarFunction& farFunction,
    const SurfaceFunction& surfaceFunction
) const
{
    return
        transportProperty
        (
            subMesh,
            &surfaceFunction::kinematicDiffusivityCoeffs,
            farFunction,
            surfaceFunction
        );
}


// ************************************************************************* //
