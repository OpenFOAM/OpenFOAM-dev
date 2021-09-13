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

#include "OSspecific.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::chemistryModel<ThermoType>::logFile(const word& name) const
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
inline Foam::label Foam::chemistryModel<ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template<class ThermoType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::chemistryModel<ThermoType>::RR()
{
    return RR_;
}


template<class ThermoType>
inline const Foam::multiComponentMixture<ThermoType>&
Foam::chemistryModel<ThermoType>::mixture() const
{
    return mixture_;
}


template<class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType>>&
Foam::chemistryModel<ThermoType>::reactions() const
{
    return reactions_;
}


template<class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::chemistryModel<ThermoType>::specieThermos() const
{
    return specieThermos_;
}


template<class ThermoType>
inline Foam::label Foam::chemistryModel<ThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class ThermoType>
inline Foam::label Foam::chemistryModel<ThermoType>::nReaction() const
{
    return nReaction_;
}


template<class ThermoType>
inline Foam::scalar Foam::chemistryModel<ThermoType>::Treact() const
{
    return Treact_;
}


template<class ThermoType>
inline Foam::scalar& Foam::chemistryModel<ThermoType>::Treact()
{
    return Treact_;
}


template<class ThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::chemistryModel<ThermoType>::RR(const label i) const
{
    return RR_[i];
}

template<class ThermoType>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::chemistryModel<ThermoType>::RR(const label i)
{
    return RR_[i];
}


template<class ThermoType>
inline const Foam::PtrList<Foam::volScalarField>&
Foam::chemistryModel<ThermoType>::Y()
{
    return this->Y_;
}


template<class ThermoType>
inline const Foam::chemistryReductionMethod<ThermoType>&
Foam::chemistryModel<ThermoType>::mechRed() const
{
    return mechRed_;
}


template<class ThermoType>
inline void Foam::chemistryModel<ThermoType>::setActive(const label i)
{
    this->thermo().composition().setActive(i);
}


template<class ThermoType>
inline bool Foam::chemistryModel<ThermoType>::active(const label i) const
{
    return this->thermo().composition().active(i);
}


template<class ThermoType>
inline void Foam::chemistryModel<ThermoType>::setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}


template<class ThermoType>
inline Foam::label Foam::chemistryModel<ThermoType>::sToc
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
inline Foam::label Foam::chemistryModel<ThermoType>::cTos
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
