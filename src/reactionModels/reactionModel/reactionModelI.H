/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

inline const Foam::fvMesh& Foam::reactionModel::mesh() const
{
    return mesh_;
}


inline const Foam::fluidMulticomponentThermo&
Foam::reactionModel::thermo() const
{
    return thermo_;
}


inline const Foam::compressibleMomentumTransportModel&
Foam::reactionModel::turbulence() const
{
    return turb_;
}


inline const Foam::volScalarField& Foam::reactionModel::rho() const
{
    return turbulence().rho();
}


inline Foam::tmp<Foam::surfaceScalarField> Foam::reactionModel::phi() const
{
    return turbulence().alphaRhoPhi();
}


inline const Foam::dictionary& Foam::reactionModel::coeffs() const
{
    return coeffs_;
}

// ************************************************************************* //
