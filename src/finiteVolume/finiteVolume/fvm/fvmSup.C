/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const DimensionedField<Type, volMesh>& su,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVolume*su.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.source() -= mesh.V()*su.primitiveField();

    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<VolField<Type>>& tsu,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


template<class Type>
Foam::zeroField
Foam::fvm::Su
(
    const zero&,
    const VolField<Type>& vf
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const volScalarField::Internal& sp,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVolume*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*sp.primitiveField();

    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField::Internal>& tsp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField>& tsp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const dimensionedScalar& sp,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVolume*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*sp.value();

    return tfvm;
}


template<class Type>
Foam::zeroField
Foam::fvm::Sp
(
    const zero&,
    const VolField<Type>&
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const volScalarField::Internal& susp,
    const VolField<Type>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVolume*susp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*max(susp.primitiveField(), scalar(0));

    fvm.source() -= mesh.V()*min(susp.primitiveField(), scalar(0))
        *vf.primitiveField();

    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField::Internal>& tsusp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField>& tsusp,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}


template<class Type>
Foam::zeroField
Foam::fvm::SuSp
(
    const zero&,
    const VolField<Type>& vf
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::S
(
    const Pair<tmp<volScalarField::Internal>>& s,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(s[1], vf) + s[0];
    return tfvm;
}


// ************************************************************************* //
