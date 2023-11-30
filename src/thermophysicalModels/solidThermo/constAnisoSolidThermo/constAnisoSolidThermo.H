/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::constAnisoSolidThermo

Description
    Uniform or non-uniform constant anisotropic solid thermodynamic properties

    Each physical property can specified as either \c uniform in which case the
    value entry is read, \c zonal in which case the value entry and zone list
    are read or \c file in which case the field file in read from the constant
    directory.  The thermal conductivity \c Kappa is anisotropic and read as a
    diagonal tensor or diagonal tensor field provided in the form of a vector
    or vector field.

Usage
    Example of uniform constant solid properties specification:
    \verbatim
        thermoType          constAnisoSolidThermo;

        rho
        {
            type        uniform;
            value       8940;
        }

        Cv
        {
            type        uniform;
            value       385;
        }

        Kappa
        {
            type        uniform;
            value       (380 100 100);
        }
    \endverbatim

    Example of zonal constant solid properties specification where Kappa is
    different in different zones:
    \verbatim
        thermoType          constSolidThermo;

        rho
        {
            type        uniform;
            value       8940;
        }

        Cv
        {
            type        uniform;
            value       385;
        }

        Kappa
        {
            type        zonal;
            value       (380 380 380);

            zones
            {
                heater      (560 560 560);
                insulation  (10  100 100);
            }
        }
    \endverbatim

    Example of non-uniform constant solid properties specification:
    \verbatim
        thermoType          constAnisoSolidThermo;

        rho
        {
            type        file;
        }

        Cv
        {
            type        file;
        }

        Kappa
        {
            type        file;
        }
    \endverbatim
    where each of the field files are read from the constant directory.

SourceFiles
    constAnisoSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constAnisoSolidThermo_H
#define constAnisoSolidThermo_H

#include "constSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class constAnisoSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class constAnisoSolidThermo
:
    public constSolidThermo
{
    // Private data

        //- Anisotropic thermal conductivity [W/m/K]
        volVectorField Kappa_;


public:

    //- Runtime type information
    TypeName("constAnisoSolidThermo");


    // Constructors

        //- Construct from mesh and phase name
        constAnisoSolidThermo
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~constAnisoSolidThermo();


    // Member Functions

        // Transport state

            //- Return false as the thermal conductivity is anisotropic
            virtual bool isotropic() const
            {
                return false;
            }

            //- Isotropic thermal conductivity [W/m/K]
            //  Not implemented
            virtual const volScalarField& kappa() const;

            //- Anisotropic thermal conductivity [W/m/K]
            virtual const volVectorField& Kappa() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
