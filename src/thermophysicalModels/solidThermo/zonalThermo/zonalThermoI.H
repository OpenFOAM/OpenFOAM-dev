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

#include "zonalThermo.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::nil Foam::zonalThermo::implementation::Yslicer() const
{
    return nil();
}


inline Foam::label Foam::zonalThermo::implementation::cellComposition
(
    const nil,
    const label celli
) const
{
    return ztzs_.cellZones()[celli];
}


inline Foam::label Foam::zonalThermo::implementation::patchFaceComposition
(
    const nil,
    const label patchi,
    const label patchFacei
) const
{
    const label facei = mesh().polyFacesBf()[patchi][patchFacei];

    return ztzs_.cellZones()[mesh().faceOwner()[facei]];
}


inline Foam::nil Foam::zonalThermo::implementation::Yslicer
(
    const fvSource& model,
    const volScalarField::Internal& source
) const
{
    return nil();
}


inline Foam::label Foam::zonalThermo::implementation::sourceCellComposition
(
    const nil,
    const label celli
) const
{
    return celli;
}


inline const Foam::labelUList& Foam::zonalThermo::implementation::Yslicer
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return cells;
}


inline Foam::label Foam::zonalThermo::implementation::sourceCellComposition
(
    const labelUList& cells,
    const label i
) const
{
    return ztzs_.cellZones()[cells[i]];
}


// ************************************************************************* //
