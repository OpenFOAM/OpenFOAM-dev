/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word&
Foam::diameterModels::velocityGroup::popBalName() const
{
    return popBalName_;
}


inline const Foam::volScalarField&
Foam::diameterModels::velocityGroup::f() const
{
    return f_;
}


inline const Foam::dimensionedScalar&
Foam::diameterModels::velocityGroup::formFactor() const
{
    return formFactor_;
}


inline const Foam::PtrList<Foam::diameterModels::sizeGroup>&
Foam::diameterModels::velocityGroup::sizeGroups() const
{
    return sizeGroups_;
}


inline const Foam::tmp<Foam::fv::convectionScheme<Foam::scalar>>&
Foam::diameterModels::velocityGroup::mvConvection() const
{
    return mvConvection_;
}


inline const Foam::volScalarField& Foam::diameterModels::velocityGroup::
dmdt() const
{
    return dmdt_;
}


inline Foam::volScalarField& Foam::diameterModels::velocityGroup::dmdt()
{
    return dmdt_;
}


// ************************************************************************* //
