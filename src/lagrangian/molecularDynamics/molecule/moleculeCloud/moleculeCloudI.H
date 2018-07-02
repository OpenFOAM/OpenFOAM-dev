/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::moleculeCloud::evaluatePair
(
    molecule& molI,
    molecule& molJ
)
{
    const pairPotentialList& pairPot = pot_.pairPotentials();

    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI.id();

    label idJ = molJ.id();

    const molecule::constantProperties& constPropI(constProps(idI));

    const molecule::constantProperties& constPropJ(constProps(idJ));

    List<label> siteIdsI = constPropI.siteIds();

    List<label> siteIdsJ = constPropJ.siteIds();

    List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();

    List<bool> electrostaticSitesI = constPropI.electrostaticSites();

    List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();

    List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

    forAll(siteIdsI, sI)
    {
        label idsI(siteIdsI[sI]);

        forAll(siteIdsJ, sJ)
        {
            label idsJ(siteIdsJ[sJ]);

            if (pairPotentialSitesI[sI] && pairPotentialSitesJ[sJ])
            {
                vector rsIsJ =
                    molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                if (pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    vector fsIsJ =
                        (rsIsJ/rsIsJMag)
                       *pairPot.force(idsI, idsJ, rsIsJMag);

                    molI.siteForces()[sI] += fsIsJ;

                    molJ.siteForces()[sJ] += -fsIsJ;

                    scalar potentialEnergy
                    (
                        pairPot.energy(idsI, idsJ, rsIsJMag)
                    );

                    molI.potentialEnergy() += 0.5*potentialEnergy;

                    molJ.potentialEnergy() += 0.5*potentialEnergy;

                    vector rIJ = molI.position() - molJ.position();

                    tensor virialContribution =
                        (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

                    molI.rf() += virialContribution;

                    molJ.rf() += virialContribution;
                }
            }

            if (electrostaticSitesI[sI] && electrostaticSitesJ[sJ])
            {
                vector rsIsJ =
                molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                if (rsIsJMagSq <= electrostatic.rCutSqr())
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    scalar chargeI = constPropI.siteCharges()[sI];

                    scalar chargeJ = constPropJ.siteCharges()[sJ];

                    vector fsIsJ =
                        (rsIsJ/rsIsJMag)
                       *chargeI*chargeJ*electrostatic.force(rsIsJMag);

                    molI.siteForces()[sI] += fsIsJ;

                    molJ.siteForces()[sJ] += -fsIsJ;

                    scalar potentialEnergy =
                        chargeI*chargeJ
                       *electrostatic.energy(rsIsJMag);

                    molI.potentialEnergy() += 0.5*potentialEnergy;

                    molJ.potentialEnergy() += 0.5*potentialEnergy;

                    vector rIJ = molI.position() - molJ.position();

                    tensor virialContribution =
                        (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

                    molI.rf() += virialContribution;

                    molJ.rf() += virialContribution;
                }
            }
        }
    }
}


inline bool Foam::moleculeCloud::evaluatePotentialLimit
(
    molecule& molI,
    molecule& molJ
) const
{
    const pairPotentialList& pairPot = pot_.pairPotentials();

    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI.id();

    label idJ = molJ.id();

    const molecule::constantProperties& constPropI(constProps(idI));

    const molecule::constantProperties& constPropJ(constProps(idJ));

    List<label> siteIdsI = constPropI.siteIds();

    List<label> siteIdsJ = constPropJ.siteIds();

    List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();

    List<bool> electrostaticSitesI = constPropI.electrostaticSites();

    List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();

    List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

    forAll(siteIdsI, sI)
    {
        label idsI(siteIdsI[sI]);

        forAll(siteIdsJ, sJ)
        {
            label idsJ(siteIdsJ[sJ]);

            if (pairPotentialSitesI[sI] && pairPotentialSitesJ[sJ])
            {
                vector rsIsJ =
                    molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                if (pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    // Guard against pairPot.energy being evaluated
                    // if rIJMag < small. A floating point exception will
                    // happen otherwise.

                    if (rsIsJMag < small)
                    {
                        WarningInFunction
                            << "Molecule site pair closer than "
                            << small
                            << ": mag separation = " << rsIsJMag
                            << ". These may have been placed on top of each"
                            << " other by a rounding error in mdInitialise in"
                            << " parallel or a block filled with molecules"
                            << " twice. Removing one of the molecules."
                            << endl;

                        return true;
                    }

                    // Guard against pairPot.energy being evaluated if rIJMag <
                    // rMin.  A tabulation lookup error will occur otherwise.

                    if (rsIsJMag < pairPot.rMin(idsI, idsJ))
                    {
                        return true;
                    }

                    if
                    (
                        mag(pairPot.energy(idsI, idsJ, rsIsJMag))
                      > pot_.potentialEnergyLimit()
                    )
                    {
                        return true;
                    };
                }
            }

            if (electrostaticSitesI[sI] && electrostaticSitesJ[sJ])
            {
                vector rsIsJ =
                    molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                if (pairPot.rCutMaxSqr(rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    // Guard against pairPot.energy being evaluated
                    // if rIJMag < small. A floating point exception will
                    // happen otherwise.

                    if (rsIsJMag < small)
                    {
                        WarningInFunction
                            << "Molecule site pair closer than "
                            << small
                            << ": mag separation = " << rsIsJMag
                            << ". These may have been placed on top of each"
                            << " other by a rounding error in mdInitialise in"
                            << " parallel or a block filled with molecules"
                            << " twice. Removing one of the molecules."
                            << endl;

                        return true;
                    }

                    if (rsIsJMag < electrostatic.rMin())
                    {
                        return true;
                    }

                    scalar chargeI = constPropI.siteCharges()[sI];

                    scalar chargeJ = constPropJ.siteCharges()[sJ];

                    if
                    (
                        mag(chargeI*chargeJ*electrostatic.energy(rsIsJMag))
                      > pot_.potentialEnergyLimit()
                    )
                    {
                        return true;
                    };
                }
            }
        }
    }

    return false;
}


inline Foam::vector Foam::moleculeCloud::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.sampleNormal<vector>();
}


inline Foam::vector Foam::moleculeCloud::equipartitionAngularMomentum
(
    scalar temperature,
    const molecule::constantProperties& cP
)
{
    scalar sqrtKbT = sqrt(physicoChemical::k.value()*temperature);

    if (cP.linearMolecule())
    {
        return sqrtKbT*vector
        (
            0.0,
            sqrt(cP.momentOfInertia().yy())*rndGen_.scalarNormal(),
            sqrt(cP.momentOfInertia().zz())*rndGen_.scalarNormal()
        );
    }
    else
    {
        return sqrtKbT*vector
        (
            sqrt(cP.momentOfInertia().xx())*rndGen_.scalarNormal(),
            sqrt(cP.momentOfInertia().yy())*rndGen_.scalarNormal(),
            sqrt(cP.momentOfInertia().zz())*rndGen_.scalarNormal()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::moleculeCloud::mesh() const
{
    return mesh_;
}


inline const Foam::potential& Foam::moleculeCloud::pot() const
{
    return pot_;
}


inline const Foam::List<Foam::DynamicList<Foam::molecule*>>&
    Foam::moleculeCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


inline const Foam::InteractionLists<Foam::molecule>&
    Foam::moleculeCloud::il() const
{
    return il_;
}


inline const Foam::List<Foam::molecule::constantProperties>
    Foam::moleculeCloud::constProps() const
{
    return constPropList_;
}


inline const Foam::molecule::constantProperties&
    Foam::moleculeCloud::constProps(label id) const
{
    return constPropList_[id];
}


inline Foam::Random& Foam::moleculeCloud::rndGen()
{
    return rndGen_;
}


// ************************************************************************* //
