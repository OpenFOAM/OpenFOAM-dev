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

#include "STARCDCoordinateRotation.H"
#include "unitConversions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(STARCDCoordinateRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        STARCDCoordinateRotation,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        STARCDCoordinateRotation,
        points
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::STARCDCoordinateRotation::calcTransform
(
    const scalar z,
    const scalar x,
    const scalar y
)
{
    R_ =
    (
        tensor
        (
            cos(y)*cos(z) - sin(x)*sin(y)*sin(z),
            -cos(x)*sin(z),
            sin(x)*cos(y)*sin(z) + sin(y)*cos(z),

            cos(y)*sin(z) + sin(x)*sin(y)*cos(z),
            cos(x)*cos(z),
            sin(y)*sin(z) - sin(x)*cos(y)*cos(z),

            -cos(x)*sin(y),
            sin(x),
            cos(x)*cos(y)
        )
    );

    Rtr_ = R_.T();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const dictionary& dict
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    const vector rotation(dict.lookup<vector>("rotation", unitDegrees));

    if (dict.found("degrees"))
    {
        FatalIOErrorInFunction(dict)
            << "Angle units are no longer specified with a 'degrees' entry. "
            << "Instead, the 'rotation' can have its units specified directly."
            << " e.g., rotation (30 45 60) [deg]."
            << exit(FatalIOError);
    }

    calcTransform
    (
        rotation.component(vector::X),
        rotation.component(vector::Y),
        rotation.component(vector::Z)
    );
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const dictionary& dict,
    const UList<vector>& points
)
:
    STARCDCoordinateRotation(dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::STARCDCoordinateRotation::transform(const vector& v) const
{
    return (R_ & v);
}


Foam::vector Foam::STARCDCoordinateRotation::invTransform(const vector& v) const
{
    return (Rtr_ & v);
}


Foam::tmp<Foam::vectorField> Foam::STARCDCoordinateRotation::transform
(
    const vectorField& vf
) const
{
    NotImplemented;
    return tmp<vectorField>(nullptr);
}


Foam::tmp<Foam::vectorField> Foam::STARCDCoordinateRotation::invTransform
(
    const vectorField& vf
) const
{
    NotImplemented;
    return tmp<vectorField>(nullptr);
}


Foam::tmp<Foam::tensorField> Foam::STARCDCoordinateRotation::transform
(
    const tensorField& tf
) const
{
    NotImplemented;
    return tmp<tensorField>(nullptr);
}


Foam::tensor Foam::STARCDCoordinateRotation::transform
(
    const vector& p,
    const tensor& t
) const
{
    return (R_ & t & Rtr_);
}


Foam::tmp<Foam::symmTensorField> Foam::STARCDCoordinateRotation::
transformDiagTensor
(
    const vectorField& vf
) const
{
    tmp<symmTensorField> tfld(new symmTensorField(vf.size()));
    symmTensorField& fld = tfld.ref();

    forAll(fld, i)
    {
        fld[i] = transformVectorDiagTensor(R_, vf[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::STARCDCoordinateRotation::transformDiagTensor
(
    const vector& p,
    const vector& v
) const
{
    return transformVectorDiagTensor(R_, v);
}


void Foam::STARCDCoordinateRotation::write(Ostream& os) const
{
     writeEntry(os, "e1", e1());
     writeEntry(os, "e2", e2());
     writeEntry(os, "e3", e3());
}


// ************************************************************************* //
