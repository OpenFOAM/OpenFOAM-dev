/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "URhoFluidMulticomponentThermo.H"
#include "uRhoMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::URhoFluidMulticomponentThermo<BaseThermo>::URhoFluidMulticomponentThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    RhoFluidThermo<BaseThermo>(mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::URhoFluidMulticomponentThermo<BaseThermo>::
~URhoFluidMulticomponentThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::URhoFluidMulticomponentThermo<BaseThermo>::Phi() const
{
    return this->volScalarFieldMixtureProperty
    (
        "Phi",
        dimless,
        &BaseThermo::mixtureType::Phi
    );
}


template<class BaseThermo>
Foam::PtrList<Foam::volScalarField::Internal>
Foam::URhoFluidMulticomponentThermo<BaseThermo>::prompt() const
{
    return BaseThermo::mixtureType::prompt(this->Y());
}


template<class BaseThermo>
void Foam::URhoFluidMulticomponentThermo<BaseThermo>::reset
(
    volScalarField& b,
    volScalarField& c,
    const PtrList<volScalarField>& Yb,
    const volScalarField& heb
)
{
    BaseThermo::mixtureType::reset(b, this->Y(), c, Yb);

    volScalarField& heu = this->he();

    for (label i=0; i<=heu.nOldTimes(); i++)
    {
        heu.oldTimeRef(i) =
            b.oldTime(i)*heu.oldTime(i) + c.oldTime(i)*heb.oldTime(i);
    }

    this->correct();

    for (label i=0; i<=b.nOldTimes(); i++)
    {
        b.oldTimeRef(i) = 1.0;
    }

    for (label i=0; i<=c.nOldTimes(); i++)
    {
        c.oldTimeRef(i) = 1.0 - b.oldTime(i);
    }
}


// ************************************************************************* //
