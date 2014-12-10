/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "mathematicalConstants.H"
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
        objectRegistry
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::STARCDCoordinateRotation::transform(const vector& st) const
{
    return (R_ & st);
}


Foam::vector Foam::STARCDCoordinateRotation::invTransform
(
    const vector& st
) const
{
    return (Rtr_ & st);
}


Foam::tmp<Foam::vectorField> Foam::STARCDCoordinateRotation::transform
(
    const vectorField& st
) const
{
    notImplemented
    (
        "tmp<vectorField> Foam::STARCDCoordinateRotation:: "
        "transform(const vectorField& st) const"
    );
    return tmp<vectorField>(NULL);
}


Foam::tmp<Foam::vectorField> Foam::STARCDCoordinateRotation::invTransform
(
    const vectorField& st
) const
{
    notImplemented
    (
        "tmp<vectorField>  Foam::STARCDCoordinateRotation::"
        "invTransform(const vectorField& st) const"
    );
    return tmp<vectorField>(NULL);
}


const Foam::tensorField& Foam::STARCDCoordinateRotation::Tr() const
{
    notImplemented
    (
        "const tensorField& STARCDCoordinateRotatio::Tr() const"
    );
     return *reinterpret_cast<const tensorField*>(0);
}


Foam::tmp<Foam::tensorField> Foam::STARCDCoordinateRotation::transformTensor
(
    const tensorField& st
) const
{
     notImplemented
    (
        "tmp<Foam::tensorField> STARCDCoordinateRotation::transformTensor()"
    );
    return tmp<tensorField>(NULL);
}


Foam::tensor Foam::STARCDCoordinateRotation::transformTensor
(
    const tensor& st
) const
{
    return (R_ & st & Rtr_);
}


Foam::tmp<Foam::tensorField> Foam::STARCDCoordinateRotation::transformTensor
(
    const tensorField& st,
    const labelList& cellMap
) const
{
    notImplemented
    (
        "tmp<Foam::tensorField> STARCDCoordinateRotation::transformTensor "
        " const tensorField& st,"
        " const labelList& cellMap "
        ") const"
    );
    return tmp<tensorField>(NULL);
}


Foam::tmp<Foam::symmTensorField> Foam::STARCDCoordinateRotation::
transformVector
(
    const vectorField& st
) const
{
    tmp<symmTensorField> tfld(new symmTensorField(st.size()));
    symmTensorField& fld = tfld();

    forAll(fld, i)
    {
        fld[i] = transformPrincipal(R_, st[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::STARCDCoordinateRotation::transformVector
(
    const vector& st
) const
{
    return transformPrincipal(R_, st);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::STARCDCoordinateRotation::calcTransform
(
    const scalar rotZ,
    const scalar rotX,
    const scalar rotY,
    const bool inDegrees
)
{
    scalar x = rotX;
    scalar y = rotY;
    scalar z = rotZ;

    if (inDegrees)
    {
        x *= constant::mathematical::pi/180.0;
        y *= constant::mathematical::pi/180.0;
        z *= constant::mathematical::pi/180.0;
    }

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

Foam::STARCDCoordinateRotation::STARCDCoordinateRotation()
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const vector& rotZrotXrotY,
    const bool inDegrees
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    calcTransform
    (
        rotZrotXrotY.component(vector::X),
        rotZrotXrotY.component(vector::Y),
        rotZrotXrotY.component(vector::Z),
        inDegrees
    );
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const scalar rotZ,
    const scalar rotX,
    const scalar rotY,
    const bool inDegrees
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    calcTransform(rotZ, rotX, rotY, inDegrees);
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const dictionary& dict
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    vector rotation(dict.lookup("rotation"));

    calcTransform
    (
        rotation.component(vector::X),
        rotation.component(vector::Y),
        rotation.component(vector::Z),
        dict.lookupOrDefault("degrees", true)
    );
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const dictionary& dict,
    const objectRegistry&
)
{
    vector rotation(dict.lookup("rotation"));

    calcTransform
    (
        rotation.component(vector::X),
        rotation.component(vector::Y),
        rotation.component(vector::Z),
        dict.lookupOrDefault("degrees", true)
    );
}


void Foam::STARCDCoordinateRotation::write(Ostream& os) const
{
     os.writeKeyword("e1") << e1() << token::END_STATEMENT << nl;
     os.writeKeyword("e2") << e2() << token::END_STATEMENT << nl;
     os.writeKeyword("e3") << e3() << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
