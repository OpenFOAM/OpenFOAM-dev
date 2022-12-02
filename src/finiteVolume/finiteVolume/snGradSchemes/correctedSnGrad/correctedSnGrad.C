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

\*---------------------------------------------------------------------------*/

#include "correctedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::correctedSnGrad<Type>::~correctedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fv::correctedSnGrad<Type>::fullGradCorrection
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct SurfaceField<Type>
    tmp<SurfaceField<Type>> tssf =
        linear<typename outerProduct<vector, Type>::type>(mesh).dotInterpolate
        (
            mesh.nonOrthCorrectionVectors(),
            gradScheme<Type>::New
            (
                mesh,
                mesh.schemes().grad("grad(" + vf.name() + ')')
            )().grad(vf, "grad(" + vf.name() + ')')
        );
    tssf.ref().rename("snGradCorr(" + vf.name() + ')');

    return tssf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fv::correctedSnGrad<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct SurfaceField<Type>
    tmp<SurfaceField<Type>> tssf
    (
        SurfaceField<Type>::New
        (
            "snGradCorr("+vf.name()+')',
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    SurfaceField<Type>& ssf = tssf.ref();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        ssf.replace
        (
            cmpt,
            correctedSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// ************************************************************************* //
