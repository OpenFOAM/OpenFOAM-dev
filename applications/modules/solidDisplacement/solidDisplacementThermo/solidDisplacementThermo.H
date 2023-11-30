/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::solidDisplacementThermo

Description
    Fundamental solid thermodynamic properties

SourceFiles
    solidDisplacementThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidDisplacementThermo_H
#define solidDisplacementThermo_H

#include "constSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidDisplacementThermo Declaration
\*---------------------------------------------------------------------------*/

class solidDisplacementThermo
:
    public constSolidThermo
{
    // Private data

        //- Switch to enable plane stress
        Switch planeStress_;

        //- Switch to enable thermal stress
        Switch thermalStress_;

        //- Youngs modulus [Pa]
        volScalarField E_;

        //- Poisson's ratio []
        volScalarField nu_;

        //- Volumetric thermal expansion coefficient [1/T]
        volScalarField alphav_;


public:

    //- Runtime type information
    TypeName("solidDisplacementThermo");


    // Constructors

        //- Construct from mesh and phase name
        solidDisplacementThermo
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidDisplacementThermo();


    // Member Functions

        //- Returns true to enable plane stress
        bool planeStress() const
        {
            return planeStress_;
        }

        //- Returns true to enable thermal stress
        bool thermalStress() const
        {
            return thermalStress_;
        }


        // Mechanical state

            //- Youngs modulus [Pa]
            virtual const volScalarField& E() const;

            //- Youngs modulus for a patch [Pa]
            virtual const scalarField& E(const label patchi) const;

            //- Poisson's ratio []
            virtual const volScalarField& nu() const;

            //- Poisson's ratio for a patch[]
            virtual const scalarField& nu(const label patchi) const;

            //- Volumetric thermal expansion coefficient [1/T]
            virtual const volScalarField& alphav() const;

            //- Volumetric thermal expansion coefficient for a patch [1/T]
            virtual const scalarField& alphav(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
