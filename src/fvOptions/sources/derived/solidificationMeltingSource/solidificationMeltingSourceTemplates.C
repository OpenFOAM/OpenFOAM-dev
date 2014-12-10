/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "fvMatrices.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::solidificationMeltingSource::apply
(
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn
)
{
    if (!solveField(eqn.psi().name()))
    {
        return;
    }

    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField Cp(this->Cp());

    update(Cp);

    dimensionedScalar L("L", dimEnergy/dimMass, L_);

    // contributions added to rhs of solver equation
    if (eqn.psi().dimensions() == dimTemperature)
    {
        // isothermal phase change - only include time derivative
//        eqn -= L/Cp*(fvc::ddt(rho, alpha1_) + fvc::div(phi, alpha1_));
        eqn -= L/Cp*(fvc::ddt(rho, alpha1_));
    }
    else
    {
        // isothermal phase change - only include time derivative
//        eqn -= L*(fvc::ddt(rho, alpha1_) + fvc::div(phi, alpha1_));
        eqn -= L*(fvc::ddt(rho, alpha1_));
    }
}


// ************************************************************************* //
