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

#include "MPLICU.H"
#include "slicedSurfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MPLICU, 0);

    surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<MPLICU>
        addMPLICUScalarMeshFluxConstructorToTable_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::MPLICU::interpolate
(
    const VolField<scalar>& vf
) const
{
    tmp<surfaceScalarField> tvff(upwind<scalar>(mesh(), phi_).interpolate(vf));

    scalarField splicedTvff
    (
        slicedSurfaceScalarField
        (
            IOobject
            (
                "splicedTvff",
                mesh().time().name(),
                mesh()
            ),
            tvff,
            false
        ).splice()
    );

    return surfaceAlpha(vf, phi_, splicedTvff, false, 1e-6);
}

// ************************************************************************* //
