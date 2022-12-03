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

#include "fvcReconstruct.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSurfaceIntegrate.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<typename outerProduct<vector, Type>::type>>
reconstruct
(
    const SurfaceField<Type>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());

    tmp<VolField<GradType>> treconField
    (
        VolField<GradType>::New
        (
            "volIntegrate("+ssf.name()+')',
            mesh,
            dimensioned<GradType>("0", ssf.dimensions()/dimArea, Zero),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );

    if (!mesh.nGeometricD())
    {
        return treconField;
    }

    treconField.ref() = inv(surfaceSum(SfHat*mesh.Sf()))&surfaceSum(SfHat*ssf),

    treconField.ref().correctBoundaryConditions();

    return treconField;
}


template<class Type>
tmp<VolField<typename outerProduct<vector, Type>::type>>
reconstruct
(
    const tmp<SurfaceField<Type>>& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<VolField<GradType>> tvf
    (
        fvc::reconstruct(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
