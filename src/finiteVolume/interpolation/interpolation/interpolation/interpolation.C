/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "interpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationBase<Type>::interpolationBase(const VolField<Type>& psi)
:
    psi_(psi),
    mesh_(psi.mesh())
{}


template<class Type>
Foam::interpolationBase<Type>::interpolationBase
(
    const interpolationBase<Type>& i
)
:
    psi_(i.psi_),
    mesh_(i.mesh_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::interpolation<Type>> Foam::interpolation<Type>::New
(
    const word& interpolationType,
    const VolField<Type>& psi
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(interpolationType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation type " << interpolationType
            << " for field " << psi.name() << nl << nl
            << "Valid interpolation types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<interpolation<Type>>(cstrIter()(psi));
}


template<class Type>
Foam::autoPtr<Foam::interpolation<Type>> Foam::interpolation<Type>::New
(
    const dictionary& interpolationSchemes,
    const VolField<Type>& psi
)
{
    return New(word(interpolationSchemes.lookup(psi.name())), psi);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationBase<Type>::~interpolationBase()
{}


template<class Type>
Foam::interpolationGradBase<Type>::~interpolationGradBase()
{}


template<class Type>
Foam::interpolation<Type>::~interpolation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolationBase<Type>::interpolate
(
    const barycentric& coordinates,
    const tetIndices& tetIs,
    const label facei
) const
{
    return
        interpolate
        (
            tetIs.tet(this->mesh_).barycentricToPoint(coordinates),
            tetIs.cell(),
            facei
        );
}


template<class Type>
Foam::interpolationGradType<Type>
Foam::interpolationGradBase<Type>::interpolateGrad
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    FatalErrorInFunction
        << "Interpolation method " << type() << " for field "
        << this->psi_.name() << " does not support gradient interpolation"
        << exit(FatalError);

    return pTraits<interpolationGradType<Type>>::nan;
}


template<class Type>
Foam::interpolationGradType<Type>
Foam::interpolationGradBase<Type>::interpolateGrad
(
    const barycentric& coordinates,
    const tetIndices& tetIs,
    const label facei
) const
{
    return
        interpolateGrad
        (
            tetIs.tet(this->mesh_).barycentricToPoint(coordinates),
            tetIs.cell(),
            facei
        );
}


template<class Type, class InterpolationType>
Foam::tmp<Foam::Field<Type>>
Foam::fieldInterpolationBase<Type, InterpolationType>::interpolate
(
    const vectorField& position,
    const labelList& celli,
    const labelList& facei
) const
{
    tmp<Field<Type>> tField(new Field<Type>(position.size()));

    Field<Type>& field = tField.ref();

    forAll(field, i)
    {
        field[i] =
            static_cast<const InterpolationType&>(*this).interpolate
            (
                position[i],
                celli[i],
                isNull(facei) ? -1 : facei[i]
            );
    }

    return tField;
}


template<class Type, class InterpolationType>
Foam::tmp<Foam::Field<Type>>
Foam::fieldInterpolationBase<Type, InterpolationType>::interpolate
(
    const Field<barycentric>& coordinates,
    const labelList& celli,
    const labelList& tetFacei,
    const labelList& tetPti,
    const labelList& facei
) const
{
    tmp<Field<Type>> tField(new Field<Type>(coordinates.size()));

    Field<Type>& field = tField.ref();

    forAll(field, i)
    {
        field[i] =
            static_cast<const InterpolationType&>(*this).interpolate
            (
                coordinates[i],
                tetIndices(celli[i], tetFacei[i], tetPti[i]),
                isNull(facei) ? -1 : facei[i]
            );
    }

    return tField;
}


template<class Type, class InterpolationType>
Foam::tmp<Foam::Field<Foam::interpolationGradType<Type>>>
Foam::fieldInterpolationGradBase<Type, InterpolationType>::interpolateGrad
(
    const vectorField& position,
    const labelList& celli,
    const labelList& facei
) const
{
    tmp<Field<interpolationGradType<Type>>> tField
    (
        new Field<interpolationGradType<Type>>(position.size())
    );

    Field<interpolationGradType<Type>>& field = tField.ref();

    forAll(field, i)
    {
        field[i] =
            static_cast<const InterpolationType&>(*this).interpolateGrad
            (
                position[i],
                celli[i],
                isNull(facei) ? -1 : facei[i]
            );
    }

    return tField;
}


template<class Type, class InterpolationType>
Foam::tmp<Foam::Field<Foam::interpolationGradType<Type>>>
Foam::fieldInterpolationGradBase<Type, InterpolationType>::interpolateGrad
(
    const Field<barycentric>& coordinates,
    const labelList& celli,
    const labelList& tetFacei,
    const labelList& tetPti,
    const labelList& facei
) const
{
    tmp<Field<interpolationGradType<Type>>> tField
    (
        new Field<interpolationGradType<Type>>(coordinates.size())
    );

    Field<interpolationGradType<Type>>& field = tField.ref();

    forAll(field, i)
    {
        field[i] =
            static_cast<const InterpolationType&>(*this).interpolateGrad
            (
                coordinates[i],
                tetIndices(celli[i], tetFacei[i], tetPti[i]),
                isNull(facei) ? -1 : facei[i]
            );
    }

    return tField;
}


// ************************************************************************* //
