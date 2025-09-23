/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "DimensionedFieldListSlicer.H"
#include "ubPsiMulticomponentThermo.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::ubPsiMulticomponentThermo::containsSpecie
(
    const word& specieName
) const
{
    return species().found(specieName);
}


inline Foam::volScalarField& Foam::ubPsiMulticomponentThermo::Y
(
    const word& specieName
)
{
    return Y()[species()[specieName]];
}


inline const Foam::volScalarField& Foam::ubPsiMulticomponentThermo::Y
(
    const word& specieName
) const
{
    return Y()[species()[specieName]];
}


inline Foam::volScalarFieldListSlicer
Foam::ubPsiMulticomponentThermo::implementation::Yslicer() const
{
    return volScalarFieldListSlicer(mesh(), Y_);
}


inline Foam::scalarFieldListSlice
Foam::ubPsiMulticomponentThermo::implementation::cellComposition
(
    const volScalarFieldListSlicer& Yslicer,
    const label celli
) const
{
    return Yslicer.slice(celli);
}


inline Foam::scalarFieldListSlice
Foam::ubPsiMulticomponentThermo::implementation::patchFaceComposition
(
    const volScalarFieldListSlicer& Yslicer,
    const label patchi,
    const label facei
) const
{
    return Yslicer.patchSlice(patchi, facei);
}


inline Foam::DimensionedFieldListAndSlicer<Foam::scalar, Foam::volMesh>
Foam::ubPsiMulticomponentThermo::implementation::Yslicer
(
    const fvSource& model,
    const volScalarField::Internal& source
) const
{
    PtrList<DimensionedField<scalar, volMesh>> result(Y_.size());

    forAll(Y_, i)
    {
        result.set
        (
            i,
            Y_[i]
           .sources()[model.name()]
           .sourceValue(model, source)
           .ptr()
        );
    }

    return result;
}


inline Foam::scalarFieldListSlice
Foam::ubPsiMulticomponentThermo::implementation::sourceCellComposition
(
    const DimensionedFieldListAndSlicer<scalar, volMesh>& Yslicer,
    const label i
) const
{
    return Yslicer.slice(i);
}


inline Foam::PtrList<Foam::scalarField>
Foam::ubPsiMulticomponentThermo::implementation::Yslicer
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    PtrList<scalarField> result(Y_.size());

    forAll(Y_, i)
    {
        result.set
        (
            i,
            Y_[i]
           .sources()[model.name()]
           .sourceValue(model, source, cells)
           .ptr()
        );
    }

    return result;
}


inline Foam::scalarFieldListSlice
Foam::ubPsiMulticomponentThermo::implementation::sourceCellComposition
(
    const PtrList<scalarField>& Yslicer,
    const label i
) const
{
    return scalarFieldListSlice(Yslicer, i);
}


// ************************************************************************* //
