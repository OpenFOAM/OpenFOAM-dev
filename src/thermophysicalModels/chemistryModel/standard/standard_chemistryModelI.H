/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

inline Foam::autoPtr<Foam::OFstream>
Foam::chemistryModels::standard::logFile(const word& name) const
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


inline Foam::label Foam::chemistryModels::standard::nSpecie() const
{
    return nSpecie_;
}


inline void Foam::chemistryModels::standard::setNSpecie(const label newNs)
{
    nSpecie_ = newNs;
}


inline Foam::label Foam::chemistryModels::standard::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::chemistryModels::standard::Y() const
{
    return Yvf_;
}


inline bool Foam::chemistryModels::standard::reduction() const
{
    return reduction_;
}


inline Foam::label Foam::chemistryModels::standard::sToc
(
    const label si
) const
{
    if (reduction_)
    {
        return sToc_[si];
    }
    else
    {
        return si;
    }
}


inline Foam::label Foam::chemistryModels::standard::cTos
(
    const label ci
) const
{
    if (reduction_)
    {
        return cTos_[ci];
    }
    else
    {
        return ci;
    }
}


// ************************************************************************* //
