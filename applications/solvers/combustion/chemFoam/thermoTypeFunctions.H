/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

template<class ThermoType>
scalarList W(const rhoReactionThermo& thermo)
{
    const PtrList<ThermoType>& specieData =
        dynamicCast<const reactingMixture<ThermoType>>(thermo)
       .speciesData();

    scalarList W(specieData.size());

    forAll(specieData, i)
    {
        W[i] = specieData[i].W();
    }

    return W;
}


template<class ThermoType>
scalar h0
(
    const rhoReactionThermo& thermo,
    const scalarList& Y,
    const scalar p,
    const scalar T
)
{
    const PtrList<ThermoType>& specieData =
        dynamic_cast<const reactingMixture<ThermoType>&>(thermo)
       .speciesData();

    scalar h0 = 0;
    forAll(Y, i)
    {
        h0 += Y[i]*specieData[i].Hs(p, T);
    }

    return h0;
}


scalarList W(const rhoReactionThermo& thermo)
{
    if (isA<reactingMixture<gasHThermoPhysics>>(thermo))
    {
        return W<gasHThermoPhysics>(thermo);
    }
    else if (isA<reactingMixture<constFluidHThermoPhysics>>(thermo))
    {
        return W<constFluidHThermoPhysics>(thermo);
    }
    else
    {
        FatalErrorInFunction
            << "Thermodynamics type " << thermo.type()
            << " not supported by chemFoam"
            << exit(FatalError);

        return scalarList::null();
    }
}


scalar h0
(
    const rhoReactionThermo& thermo,
    const scalarList& Y,
    const scalar p,
    const scalar T
)
{
    if (isA<reactingMixture<gasHThermoPhysics>>(thermo))
    {
        return h0<gasHThermoPhysics>(thermo, Y, p, T);
    }
    else if (isA<reactingMixture<constFluidHThermoPhysics>>(thermo))
    {
        return h0<constFluidHThermoPhysics>(thermo, Y, p, T);
    }
    else
    {
        FatalErrorInFunction
            << "Thermodynamics type " << thermo.type()
            << " not supported by chemFoam"
            << exit(FatalError);

        return 0;
    }
}


// ************************************************************************* //
