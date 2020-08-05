/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "EulerImplicit.H"
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::EulerImplicit
(
    const fluidReactionThermo& thermo
)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("EulerImplicitCoeffs")),
    cTauChem_(coeffsDict_.lookup<scalar>("cTauChem")),
    cTp_(this->nEqns()),
    R_(this->nEqns()),
    J_(this->nEqns()),
    E_(this->nEqns() - 2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::~EulerImplicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::EulerImplicit<ChemistryModel>::solve
(
    scalar& p,
    scalar& T,
    scalarField& c,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    const label nSpecie = this->nSpecie();

    // Map the composition, temperature and pressure into cTp
    for (int i=0; i<nSpecie; i++)
    {
        cTp_[i] = max(0, c[i]);
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie + 1] = p;

    // Calculate the reaction rate and Jacobian
    this->jacobian(0, cTp_, li, R_, J_);

    // Calculate the stable/accurate time-step
    scalar tMin = great;
    const scalar cTot = sum(c);

    for (label i=0; i<nSpecie; i++)
    {
        if (R_[i] < -small)
        {
            tMin = min(tMin, -(cTp_[i] + small)/R_[i]);
        }
        else
        {
            tMin = min
            (
                tMin,
                max(cTot - cTp_[i], 1e-5)/max(R_[i], small)
            );
        }
    }

    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);

    // Assemble the Euler implicit matrix for the composition
    scalarField& source = E_.source();
    for (label i=0; i<nSpecie; i++)
    {
        E_(i, i) = 1/deltaT - J_(i, i);
        source[i] = R_[i] + E_(i, i)*cTp_[i];

        for (label j=0; j<nSpecie; j++)
        {
            if (i != j)
            {
                E_(i, j) = -J_(i, j);
                source[i] += E_(i, j)*cTp_[j];
            }
        }
    }

    // Solve for the new composition
    scalarField::subField(cTp_, nSpecie) = E_.LUsolve();

    // Limit the composition and transfer back into c
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0, cTp_[i]);
    }

    // Euler explicit integrate the temperature.
    // Separating the integration of temperature from composition
    // is significantly more stable for exothermic systems
    T += deltaT*R_[nSpecie];
}


// ************************************************************************* //
