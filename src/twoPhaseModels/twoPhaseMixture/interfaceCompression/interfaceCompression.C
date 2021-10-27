/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "interfaceCompression.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompressionNew, 0);

    surfaceInterpolationScheme<scalar>::
        addMeshFluxConstructorToTable<interfaceCompressionNew>
        addinterfaceCompressionScalarMeshFluxConstructorToTable_;

    const wordHashSet compressionSchemes
    {
        "interfaceCompression",
        "noInterfaceCompression",
        "PLIC",
        "PLICU",
        "MPLIC",
        "MPLICU"
    };
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceCompressionNew::interpolate
(
    const volScalarField& vf
) const
{
    const surfaceScalarField& nHatf =
        mesh().lookupObject<const surfaceScalarField>
        (
            "nHatf"
        );

    const surfaceScalarField vff
    (
        linear<scalar>(mesh()).interpolate(vf)
    );

    surfaceScalarField vfc
    (
        cAlpha_*sign(phi_)*vff*(1 - vff)*nHatf/mesh().magSf()
    );

    surfaceScalarField::Boundary& vfcBf = vfc.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(vfc.boundaryField(), patchi)
    {
        fvsPatchScalarField& vfcp = vfcBf[patchi];

        if (!vfcp.coupled())
        {
            vfcp == 0;
        }
    }

    tmp<surfaceScalarField> tvff(tScheme_().interpolate(vf) + vfc);
    tvff.ref().maxMin(0, 1);

    return tvff;
}


// ************************************************************************* //
