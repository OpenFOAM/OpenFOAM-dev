/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "axesRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(axesRotation, 0);
    addToRunTimeSelectionTable(coordinateRotation, axesRotation, dictionary);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        axesRotation,
        objectRegistry
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::axesRotation::calcTransform
(
    const vector& axis1,
    const vector& axis2,
    const axisOrder& order
)
{
    vector a = axis1/mag(axis1);
    vector b = axis2;

    b = b - (b & a)*a;

    if (mag(b) < small)
    {
        FatalErrorInFunction
            << "axis1, axis2 appear co-linear: "
            << axis1 << ", " << axis2 << endl
            << abort(FatalError);
    }

    b = b/mag(b);
    vector c = a^b;

    tensor Rtr;
    switch (order)
    {
        case e1e2:
        {
            Rtr = tensor(a, b, c);
            break;
        }
        case e2e3:
        {
            Rtr = tensor(c, a, b);
            break;
        }
        case e3e1:
        {
            Rtr = tensor(b, c, a);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled axes specifictation" << endl
                << abort(FatalError);

            Rtr = Zero;
            break;
        }
    }

    // Global->local transformation
    Rtr_ = Rtr;

    // Local->global transformation
    R_ = Rtr.T();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::axesRotation::axesRotation()
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{}


Foam::axesRotation::axesRotation
(
    const vector& axis,
    const vector& dir
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    calcTransform(axis, dir, e3e1);
}


Foam::axesRotation::axesRotation
(
    const dictionary& dict
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    operator=(dict);
}


Foam::axesRotation::axesRotation
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    R_(sphericalTensor::I),
    Rtr_(R_)
{
    operator=(dict);
}


Foam::axesRotation::axesRotation(const tensor& R)
:
    R_(R),
    Rtr_(R_.T())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::tensorField& Foam::axesRotation::Tr() const
{
    NotImplemented;
    return NullObjectRef<tensorField>();
}


Foam::tmp<Foam::vectorField> Foam::axesRotation::transform
(
    const vectorField& st
) const
{
    return (R_ & st);
}


Foam::vector Foam::axesRotation::transform(const vector& st) const
{
    return (R_ & st);
}


Foam::tmp<Foam::vectorField> Foam::axesRotation::invTransform
(
    const vectorField& st
) const
{
    return (Rtr_ & st);
}


Foam::vector Foam::axesRotation::invTransform(const vector& st) const
{
    return (Rtr_ & st);
}


Foam::tmp<Foam::tensorField> Foam::axesRotation::transformTensor
(
    const tensorField& st
) const
{
    NotImplemented;
    return tmp<tensorField>(nullptr);
}


Foam::tensor Foam::axesRotation::transformTensor
(
    const tensor& st
) const
{
    return (R_ & st & Rtr_);
}


Foam::tmp<Foam::tensorField> Foam::axesRotation::transformTensor
(
    const tensorField& st,
    const labelList& cellMap
) const
{
    NotImplemented;
    return tmp<tensorField>(nullptr);
}


Foam::tmp<Foam::symmTensorField> Foam::axesRotation::transformVector
(
    const vectorField& st
) const
{
    tmp<symmTensorField> tfld(new symmTensorField(st.size()));
    symmTensorField& fld = tfld.ref();

    forAll(fld, i)
    {
        fld[i] = transformPrincipal(R_, st[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::axesRotation::transformVector
(
    const vector& st
) const
{
    return transformPrincipal(R_, st);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::axesRotation::operator=(const dictionary& dict)
{
    if (debug)
    {
        Pout<< "axesRotation::operator=(const dictionary&) : "
            << "assign from " << dict << endl;
    }

    vector axis1, axis2;
    axisOrder order(e3e1);

    if (dict.readIfPresent("e1", axis1) && dict.readIfPresent("e2", axis2))
    {
        order = e1e2;
    }
    else if (dict.readIfPresent("e2", axis1)&& dict.readIfPresent("e3", axis2))
    {
        order = e2e3;
    }
    else if (dict.readIfPresent("e3", axis1)&& dict.readIfPresent("e1", axis2))
    {
        order = e3e1;
    }
    else if (dict.found("axis") || dict.found("direction"))
    {
        // Both "axis" and "direction" are required
        // If one is missing the appropriate error message will be generated
        order = e3e1;
        dict.lookup("axis") >> axis1;
        dict.lookup("direction") >> axis2;
    }
    else
    {
        FatalErrorInFunction
            << "not entry of the type (e1, e2) or (e2, e3) or (e3, e1) "
            << "found "
            << exit(FatalError);
    }

    calcTransform(axis1, axis2, order);
}


void Foam::axesRotation::write(Ostream& os) const
{
     os.writeKeyword("e1") << e1() << token::END_STATEMENT << nl;
     os.writeKeyword("e2") << e2() << token::END_STATEMENT << nl;
     os.writeKeyword("e3") << e3() << token::END_STATEMENT << nl;
}


// ************************************************************************* //
