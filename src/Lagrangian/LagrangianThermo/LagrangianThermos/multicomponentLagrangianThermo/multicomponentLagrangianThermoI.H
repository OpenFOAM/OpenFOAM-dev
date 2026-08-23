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

#include "multicomponentLagrangianThermo.H"
#include "LagrangianSubFieldsFwd.H"
#include "LagrangianFieldSource.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::LagrangianScalarDynamicField&
Foam::multicomponentLagrangianThermo::Y(const label speciei)
{
    return Y()[speciei];
}


inline const Foam::LagrangianScalarDynamicField&
Foam::multicomponentLagrangianThermo::Y(const label speciei) const
{
    return Y()[speciei];
}


inline Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::multicomponentLagrangianThermo::Y
(
    const label speciei,
    const LagrangianSubMesh& subMesh
) const
{
    return subMesh.sub(Y()[speciei]);
}


inline Foam::LagrangianScalarDynamicField&
Foam::multicomponentLagrangianThermo::Y(const word& specieName)
{
    return Y()[species()[specieName]];
}


inline const Foam::LagrangianScalarDynamicField&
Foam::multicomponentLagrangianThermo::Y(const word& specieName) const
{
    return Y()[species()[specieName]];
}


inline Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::multicomponentLagrangianThermo::Y
(
    const word& specieName,
    const LagrangianSubMesh& subMesh
) const
{
    return subMesh.sub(Y()[species()[specieName]]);
}


inline Foam::label Foam::multicomponentLagrangianThermo::specieIndex
(
    const LagrangianScalarDynamicField& Yi
) const
{
    return species()[Yi.member()];
}


inline Foam::DimensionedFieldListSlicer<Foam::scalar>
Foam::multicomponentLagrangianThermo::implementation::Yslicer() const
{
    return DimensionedFieldListSlicer<scalar>(Y_);
}


inline Foam::DimensionedFieldListSlicer<Foam::scalar>
Foam::multicomponentLagrangianThermo::implementation::Yslicer
(
    const LagrangianSubMesh& subMesh
) const
{
    return
        DimensionedFieldListSlicer<scalar>
        (
            Y_,
            subMesh.size(),
            subMesh.start()
        );
}


inline
Foam::DimensionedFieldListAndSlicer<Foam::scalar, Foam::LagrangianSubMesh>
Foam::multicomponentLagrangianThermo::implementation::Yslicer
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    PtrList<DimensionedField<scalar, LagrangianSubMesh>> result(Y_.size());

    forAll(Y_, i)
    {
        result.set
        (
            i,
            Y_[i].sources()[injection.name()].value(injection, subMesh).ptr()
        );
    }

    return result;
}


inline Foam::scalarFieldListSlice
Foam::multicomponentLagrangianThermo::implementation::
injectionElementComposition
(
    const DimensionedFieldListAndSlicer<scalar, LagrangianSubMesh>& Yslicer,
    const label subi
) const
{
    return Yslicer.slice(subi);
}


// ************************************************************************* //
