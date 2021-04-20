/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "ode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ode<ChemistryModel>::ode(const fluidReactionThermo& thermo)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("odeCoeffs")),
    odeSolver_(ODESolver::New(*this, coeffsDict_)),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ode<ChemistryModel>::~ode()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::ode<ChemistryModel>::solve
(
    scalar& p,
    scalar& T,
    scalarField& c,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(cTp_);
    }

    const label nSpecie = this->nSpecie();

    // Copy the concentration, T and P to the total solve-vector
    for (int i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = p;

    if (debug)
    {
        scalarField dcTp(this->nEqns(), rootSmall);
        dcTp[nSpecie] = 300*rootSmall;
        dcTp[nSpecie+1] = 1e5*rootSmall;
        this->check(0, cTp_, dcTp, li);
    }

    odeSolver_->solve(0, deltaT, cTp_, li, subDeltaT);

    for (int i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, cTp_[i]);
    }
    T = cTp_[nSpecie];
    p = cTp_[nSpecie+1];
}


// ************************************************************************* //
