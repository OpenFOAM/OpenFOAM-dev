/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "phasePair.H"
#include "phaseSystem.H"
#include "surfaceTensionModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH
(
    const volScalarField& d
) const
{
    return
        mag(dispersed().rho() - continuous().rho())
       *mag(g())
       *sqr(d)
       /sigma();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePair::phasePair
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const bool ordered
)
:
    phasePairKey(phase1.name(), phase2.name(), ordered),
    phase1_(phase1),
    phase2_(phase2),
    g_(phase1.mesh().lookupObject<uniformDimensionedVectorField>("g"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePair::~phasePair()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phasePair::dispersed() const
{
    FatalErrorInFunction
        << "Requested dispersed phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


const Foam::phaseModel& Foam::phasePair::continuous() const
{
    FatalErrorInFunction
        << "Requested continuous phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


Foam::word Foam::phasePair::name() const
{
    word name2(second());
    name2[0] = toupper(name2[0]);
    return first() + "And" + name2;
}


Foam::word Foam::phasePair::otherName() const
{
    word name1(first());
    name1[0] = toupper(name1[0]);
    return second() + "And" + name1;
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::rho() const
{
    return phase1()*phase1().rho() + phase2()*phase2().rho();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::magUr() const
{
    return mag(phase1().U() - phase2().U());
}


Foam::tmp<Foam::volVectorField> Foam::phasePair::Ur() const
{
    return dispersed().U() - continuous().U();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Re() const
{
    return magUr()*dispersed().d()/continuous().nu();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Pr() const
{
    return
         continuous().nu()
        *continuous().thermo().Cpv()
        *continuous().rho()
        /continuous().kappa();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Eo() const
{
    return EoH(dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH1() const
{
    return
        EoH
        (
            dispersed().d()
           *cbrt(1 + 0.163*pow(Eo(), 0.757))
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH2() const
{
    return
        EoH
        (
            dispersed().d()
           /cbrt(E())
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::sigma() const
{
    return
        phase1().fluid().lookupSubModel<surfaceTensionModel>
        (
            phasePair(phase1(), phase2())
        ).sigma();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Mo() const
{
    return
        mag(g())
       *continuous().nu()
       *pow3
        (
            continuous().nu()
           *continuous().rho()
           /sigma()
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Ta() const
{
    return Re()*pow(Mo(), 0.23);
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::E() const
{
    FatalErrorInFunction
        << "Requested aspect ratio of the dispersed phase in an unordered pair"
        << exit(FatalError);

    return phase1();
}


// ************************************************************************* //
