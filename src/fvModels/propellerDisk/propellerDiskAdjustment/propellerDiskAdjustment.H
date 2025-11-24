/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

Class
    Foam::fv::propellerDiskAdjustment

Description
    Automatic adjustment option for the propellerDisk momentum source

    A proportional-integral (PI) controller is used to adjust the propeller
    rotation speed to achieve balance between the propeller thrust and the hull
    resistance.  A user-specified relaxation time is used to control the rate
    at which the porpeller speed is adjusted.

    Reference:
    \verbatim
        Nuutinen, M. (2019).
        Automated self-propulsion point search algorithm for ship performance
        CFD simulations.
        Sixth International Symposium on Marine Propulsors, SMP’19 Rome, Italy.
    \endverbatim

Usage
    Example usage:
    \verbatim
    diskSource
    {
        type            propellerDisk;

        cellZone        propeller;

        normal          (1 0 0);    // Normal direction of the propeller

        n               26.03;      // Rotation speed [1/s]

        dPropeller      0.203;      // Propeller diameter
        dHub            0.039179;   // Hub diameter

        // Automatic adjustment controls
        adjustment      yes;
        startTime       0;      // Start time of rotation speed adjustment
        deltaTStar      0.1;    // Relaxation time scaler [s]
        nFraction       0.1;    // Maximum fractional change of rotation speed
        sfc             0;      // Skin-friction correction [N]
        Tmin            1;      // Minimum thrust for beta [N]

        resistanceFraction  1;  // Defaults to 1
        resistanceDirection $normal;    // Defaults to $normal

        forces
        {
            type            forces;
            libs            ("libforces.so");
            patches         (pod pod_end napa napa_karki rako);
            log             on;
            writeControl    timeStep;
            writeInterval   1;
            CofR            (0 0 0);

            rho             rhoInf;
            rhoInf          998.8;
        }

        propellerCurve
        {
            type table;

            //   J     Kt     Kq
            values
            (
                (0.10 (0.3267 0.03748))
                (0.15 (0.3112 0.03629))
                (0.20 (0.2949 0.03500))
                (0.25 (0.2777 0.03361))
                (0.30 (0.2598 0.03210))
                (0.35 (0.2410 0.03047))
                (0.40 (0.2214 0.02871))
                (0.45 (0.2010 0.02682))
                (0.50 (0.1798 0.02479))
                (0.55 (0.1577 0.02261))
                (0.60 (0.1349 0.02027))
                (0.65 (0.1112 0.01777))
                (0.70 (0.0867 0.01509))
                (0.75 (0.0614 0.01224))
                (0.80 (0.0353 0.00921))
            );
        }
    }
    \endverbatim

SourceFiles
    propellerDiskAdjustment.C
    propellerDiskAdjustmentTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef propellerDiskAdjustment_H
#define propellerDiskAdjustment_H

#include "propellerDisk.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class propellerDiskAdjustment Declaration
\*---------------------------------------------------------------------------*/

class propellerDiskAdjustment
{
protected:

    // Protected Data

        //- Reference to the propellerDisk fvModel
        const propellerDisk& propellerDisk_;

        //- Function object to calculate the resistance force
        mutable autoPtr<functionObjects::forces> forces_;

        //- Start time for rotation speed adjustment
        scalar startTime_;

        //- Relaxation time scale for rotation speed adjustment
        scalar deltaTStar_;

        //- Maximum fractional change of rotation speed
        scalar nFraction_;

        //- Minimum thrust for beta [N]
        scalar Tmin_;

        //- Skin friction correction for model-scale ship [N]
        scalar sfc_;

        //- Resistance fraction for multi-propulsor configurations
        scalar resistanceFraction_;

        //- Resistance force direction
        vector resistanceDirection_;

        //- Self-propulsion corrected rotational speed [1/s]
        mutable uniformDimensionedScalarField n_;


    // Protected Member Functions

        //- Return the ship resistance (for self-propulsion correction)
        scalar resistance() const;

        //- Return the current rotation speed
        scalar n() const
        {
            return n_.value();
        }

        //- Correct the rotation speed from the current propulsion force
        void correctn(const scalar T) const;

        friend class propellerDisk;


private:

    // Private Member Functions

        //- Read the model coefficients
        void readCoeffs(const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        propellerDiskAdjustment
        (
            const propellerDisk&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        propellerDiskAdjustment(const propellerDiskAdjustment&) = delete;


    //- Destructor
    ~propellerDiskAdjustment()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const propellerDiskAdjustment&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
