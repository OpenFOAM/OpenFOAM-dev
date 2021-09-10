/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

template<class ThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::TDACChemistryModel<ThermoType>::logFile(const word& name) const
{
    mkDir(this->mesh().time().path()/"TDAC"/this->group());
    return autoPtr<OFstream>
    (
        new OFstream
        (
            this->mesh().time().path()/"TDAC"/this->group()/name
        )
    );
}


template<class ThermoType>
inline const Foam::PtrList<Foam::volScalarField>&
Foam::TDACChemistryModel<ThermoType>::Y()
{
    return this->Y_;
}


template<class ThermoType>
inline const Foam::chemistryReductionMethod<ThermoType>&
Foam::TDACChemistryModel<ThermoType>::mechRed() const
{
    return mechRed_;
}


template<class ThermoType>
inline void Foam::TDACChemistryModel<ThermoType>::setActive(const label i)
{
    this->thermo().composition().setActive(i);
}


template<class ThermoType>
inline bool Foam::TDACChemistryModel<ThermoType>::active(const label i) const
{
    return this->thermo().composition().active(i);
}


template<class ThermoType>
inline void Foam::TDACChemistryModel<ThermoType>::setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}


template<class ThermoType>
inline Foam::label Foam::TDACChemistryModel<ThermoType>::sToc
(
    const label si
) const
{
    if (mechRedActive_)
    {
        return sToc_[si];
    }
    else
    {
        return si;
    }
}


template<class ThermoType>
inline Foam::label Foam::TDACChemistryModel<ThermoType>::cTos
(
    const label ci
) const
{
    if (mechRedActive_)
    {
        return cTos_[ci];
    }
    else
    {
        return ci;
    }
}


// ************************************************************************* //
