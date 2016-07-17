/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

template<class CompType, class ThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::TDACChemistryModel<CompType, ThermoType>::logFile(const word& name) const
{
    mkDir(this->mesh().time().path()/"TDAC");
    return autoPtr<OFstream>
    (
        new OFstream
        (
            this->mesh().time().path()/"TDAC"/name
        )
    );
}


template<class CompType, class ThermoType>
inline Foam::PtrList<Foam::volScalarField>&
Foam::TDACChemistryModel<CompType, ThermoType>::Y()
{
    return this->Y_;
}


template<class CompType, class ThermoType>
inline Foam::autoPtr<Foam::chemistryReductionMethod<CompType, ThermoType>>&
Foam::TDACChemistryModel<CompType, ThermoType>::mechRed()
{
    return mechRed_;
}


template<class CompType, class ThermoType>
inline void Foam::TDACChemistryModel<CompType, ThermoType>::setActive
(
    const label i
)
{
    this->thermo().composition().setActive(i);
}


template<class CompType, class ThermoType>
inline bool Foam::TDACChemistryModel<CompType, ThermoType>::active
(
    const label i
) const
{
    return this->thermo().composition().active(i);
}


template<class CompType, class ThermoType>
inline void
Foam::TDACChemistryModel<CompType, ThermoType>::setNsDAC(const label newNsDAC)
{
    NsDAC_ = newNsDAC;
}


template<class CompType, class ThermoType>
inline void
Foam::TDACChemistryModel<CompType, ThermoType>::setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}


template<class CompType, class ThermoType>
inline Foam::DynamicList<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::simplifiedToCompleteIndex()
{
    return simplifiedToCompleteIndex_;
}


template<class CompType, class ThermoType>
inline Foam::Field<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::completeToSimplifiedIndex()
{
    return completeToSimplifiedIndex_;
}


template<class CompType, class ThermoType>
inline const Foam::Field<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::
completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


template<class CompType, class ThermoType>
inline Foam::Field<bool>&
Foam::TDACChemistryModel<CompType, ThermoType>::reactionsDisabled()
{
    return reactionsDisabled_;
}


template<class CompType, class ThermoType>
inline Foam::scalarField&
Foam::TDACChemistryModel<CompType, ThermoType>::completeC()
{
    return completeC_;
}


template<class CompType, class ThermoType>
inline Foam::scalarField&
Foam::TDACChemistryModel<CompType, ThermoType>::simplifiedC()
{
    return simplifiedC_;
}


template<class CompType, class ThermoType>
inline Foam::List<Foam::List<Foam::specieElement>>&
Foam::TDACChemistryModel<CompType, ThermoType>::specieComp()
{
    return specieComp_;
}


// ************************************************************************* //
