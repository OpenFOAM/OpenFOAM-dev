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
inline Foam::label Foam::TDACChemistryModel<ThermoType>::timeSteps() const
{
    return timeSteps_;
}


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
inline Foam::autoPtr<Foam::chemistryReductionMethod<ThermoType>>&
Foam::TDACChemistryModel<ThermoType>::mechRed()
{
    return mechRed_;
}


template<class ThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::TDACChemistryModel<ThermoType>::tabulationResults() const
{
    return tabulationResults_;
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
inline void Foam::TDACChemistryModel<ThermoType>::setNsDAC(const label newNsDAC)
{
    NsDAC_ = newNsDAC;
}


template<class ThermoType>
inline void Foam::TDACChemistryModel<ThermoType>::setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}


template<class ThermoType>
inline Foam::DynamicList<Foam::label>&
Foam::TDACChemistryModel<ThermoType>::simplifiedToCompleteIndex()
{
    return simplifiedToCompleteIndex_;
}


template<class ThermoType>
inline Foam::Field<Foam::label>&
Foam::TDACChemistryModel<ThermoType>::completeToSimplifiedIndex()
{
    return completeToSimplifiedIndex_;
}


template<class ThermoType>
inline const Foam::Field<Foam::label>&
Foam::TDACChemistryModel<ThermoType>::completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


template<class ThermoType>
inline Foam::Field<bool>&
Foam::TDACChemistryModel<ThermoType>::reactionsDisabled()
{
    return reactionsDisabled_;
}


template<class ThermoType>
inline Foam::scalarField& Foam::TDACChemistryModel<ThermoType>::completeC()
{
    return completeC_;
}


template<class ThermoType>
inline Foam::scalarField& Foam::TDACChemistryModel<ThermoType>::simplifiedC()
{
    return simplifiedC_;
}


template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::resetTabulationResults()
{
    forAll(tabulationResults_, tabi)
    {
        tabulationResults_[tabi] = 2;
    }
}


// ************************************************************************* //
