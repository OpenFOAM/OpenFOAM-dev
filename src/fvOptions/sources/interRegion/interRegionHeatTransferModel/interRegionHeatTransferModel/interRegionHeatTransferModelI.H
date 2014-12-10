/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::volScalarField&
Foam::fv::interRegionHeatTransferModel::htc() const
{
    return  htc_;
}


inline const Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel() const
{
    if (nbrModel_ == NULL)
    {
        FatalErrorIn("const interRegionHeatTransferModel& nbrModel() const")
            << "Neighbour model not set"
            << abort(FatalError);
    }

    return *nbrModel_;
}


inline Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel()
{
    if (nbrModel_ == NULL)
    {
        FatalErrorIn("interRegionHeatTransferModel& nbrModel()")
            << "Neighbour model not set"
            << abort(FatalError);
    }

    return *nbrModel_;
}


// ************************************************************************* //
