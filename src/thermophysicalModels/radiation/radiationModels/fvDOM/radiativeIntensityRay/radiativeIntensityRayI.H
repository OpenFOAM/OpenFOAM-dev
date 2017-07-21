/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRay::I() const
{
    return I_;
}


inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRay::qr() const
{
    return qr_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRay::qr()
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiation::
radiativeIntensityRay::qin() const
{
    return qin_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRay::qin()
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiation::
radiativeIntensityRay::qem() const
{
    return qem_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRay::qem()
{
    return qem_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRay::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRay::dAve() const
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::omega() const
{
    return omega_;
}


inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRay::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
