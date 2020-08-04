/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "fvcDdt.H"
#include "twoPhaseMixtureThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::VoFSolidificationMeltingSource::apply
(
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    update();

    const twoPhaseMixtureThermo& thermo
    (
        mesh_.lookupObject<twoPhaseMixtureThermo>
        (
            twoPhaseMixtureThermo::dictName
        )
    );

    const volScalarField CpVoF(thermo.thermo1().Cp());

    if (eqn.psi().dimensions() == dimTemperature)
    {
        eqn += L_/CpVoF*(fvc::ddt(rho, alphaSolid_));
    }
    else
    {
        eqn += L_*(fvc::ddt(rho, alphaSolid_));
    }
}


// ************************************************************************* //
