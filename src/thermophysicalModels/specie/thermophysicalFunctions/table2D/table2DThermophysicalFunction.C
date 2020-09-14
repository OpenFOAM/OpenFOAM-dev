/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "table2DThermophysicalFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{
    defineTypeNameAndDebug(table2D, 0);
    addToRunTimeSelectionTable(thermophysicalFunction, table2D, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalFunctions::table2D::table2D
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    pLow_(dict.lookup<scalar>("pLow")),
    pHigh_(dict.lookup<scalar>("pHigh")),
    Tlow_(dict.lookup<scalar>("Tlow")),
    Thigh_(dict.lookup<scalar>("Thigh")),
    values_(dict.lookup(name))
{
    if (values_.m() < 2 || values_.n() < 2)
    {
        FatalErrorInFunction
            << "Table " << nl
            << "    " << name_ << nl
            << "    has less than 2 entries in one or both dimensions."
            << exit(FatalError);
    }
    else
    {
        deltap_ = (pHigh_ - pLow_)/(values_.m() - 1);
        deltaT_ = (Thigh_ - Tlow_)/(values_.n() - 1);
    }
}


Foam::thermophysicalFunctions::table2D::table2D(const dictionary& dict)
:
    table2D("values", dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::thermophysicalFunctions::table2D::checkRange
(
    scalar p,
    scalar ndp,
    label ip,
    scalar T,
    scalar ndT,
    label iT
) const
{
    if (ndp < 0 || ip > values_.m() - 2)
    {
        FatalErrorInFunction
            << "Pressure " << p << " out of range "
            << pLow_ << " to " << pHigh_ << nl
            << "    of table " << name_
            << exit(FatalError);
    }

    if (ndT < 0 || iT > values_.n() - 2)
    {
        FatalErrorInFunction
            << "Temperature " << T << " out of range "
            << Tlow_ << " to " << Thigh_ << nl
            << "    of table " << name_
            << exit(FatalError);
    }
}


Foam::scalar Foam::thermophysicalFunctions::table2D::f(scalar p, scalar T) const
{
    const scalar ndp = (p - pLow_)/deltap_;
    const label ip = ndp;

    const scalar ndT = (T - Tlow_)/deltaT_;
    const label iT = ndT;

    checkRange(p, ndp, ip, T, ndT, iT);

    const scalar pi = pLow_ + ip*deltap_;
    const scalar lambdap = (p - pi)/deltap_;

    // Interpolate the values at Ti wrt p
    const scalar fpi =
        values_(ip, iT)
      + lambdap*(values_(ip + 1, iT) - values_(ip, iT));

    // Interpolate the values at Ti+1 wrt p
    const scalar fpip1 =
        values_(ip, iT + 1)
      + lambdap*(values_(ip + 1, iT + 1) - values_(ip, iT + 1));

    const scalar Ti = Tlow_ + iT*deltaT_;
    const scalar lambdaT = (T - Ti)/deltaT_;

    // Interpolate wrt T
    return fpi + lambdaT*(fpip1 - fpi);
}


Foam::scalar Foam::thermophysicalFunctions::table2D::
dfdp
(
    scalar p,
    scalar T
) const
{
    const scalar ndp = (p - pLow_)/deltap_;
    const label ip = ndp;

    const scalar ndT = (T - Tlow_)/deltaT_;
    const label iT = ndT;

    checkRange(p, ndp, ip, T, ndT, iT);

    const scalar dfdpi =
        (values_(ip + 1, iT) - values_(ip, iT))/deltap_;
    const scalar dfdpip1 =
        (values_(ip + 1, iT + 1) - values_(ip, iT + 1))/deltap_;

    const scalar Ti = Tlow_ + iT*deltaT_;
    const scalar lambdaT = (T - Ti)/deltaT_;

    // Interpolate wrt T
    return dfdpi + lambdaT*(dfdpip1 - dfdpi);
}


Foam::scalar Foam::thermophysicalFunctions::table2D::
dfdT
(
    scalar p,
    scalar T
) const
{
    const scalar ndp = (p - pLow_)/deltap_;
    const label ip = ndp;

    const scalar ndT = (T - Tlow_)/deltaT_;
    const label iT = ndT;

    checkRange(p, ndp, ip, T, ndT, iT);

    const scalar dfdTi =
        (values_(ip, iT + 1) - values_(ip, iT))/deltaT_;
    const scalar dfdTip1 =
        (values_(ip + 1, iT + 1) - values_(ip + 1, iT))/deltaT_;

    const scalar pi = pLow_ + ip*deltap_;
    const scalar lambdap = (p - pi)/deltap_;

    // Interpolate wrt p
    return dfdTi + lambdap*(dfdTip1 - dfdTi);
}


void Foam::thermophysicalFunctions::table2D::write(Ostream& os) const
{
    writeEntry(os, "pLow", Tlow_);
    writeEntry(os, "pHigh", Thigh_);
    writeEntry(os, "Tlow", Tlow_);
    writeEntry(os, "Thigh", Thigh_);
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
