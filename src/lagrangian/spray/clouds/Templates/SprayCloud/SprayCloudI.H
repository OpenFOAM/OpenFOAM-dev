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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::SprayCloud<CloudType>&
Foam::SprayCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::AtomizationModel<Foam::SprayCloud<CloudType> >&
Foam::SprayCloud<CloudType>::atomization() const
{
    return atomizationModel_;
}


template<class CloudType>
inline Foam::AtomizationModel<Foam::SprayCloud<CloudType> >&
Foam::SprayCloud<CloudType>::atomization()
{
    return atomizationModel_();
}


template<class CloudType>
inline const Foam::BreakupModel<Foam::SprayCloud<CloudType> >&
Foam::SprayCloud<CloudType>::breakup() const
{
    return breakupModel_;
}


template<class CloudType>
inline Foam::BreakupModel<Foam::SprayCloud<CloudType> >&
Foam::SprayCloud<CloudType>::breakup()
{
    return breakupModel_();
}


template<class CloudType>
inline Foam::scalar Foam::SprayCloud<CloudType>::averageParcelMass() const
{
    return averageParcelMass_;
}


// ************************************************************************* //
