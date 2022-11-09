/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::interfaceProperties

Description
    Contains the interface properties. Surface tension, interface curvature,
    and functionality to correct alpha contact angle boundary conditions.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private Data

        //- Keep a reference to the phaseProperties dictionary
        const dictionary& phasePropertiesDict_;

        //- Volume fraction of phase 1
        volScalarField& alpha1_;

        //- Volume fraction of phase 2
        volScalarField& alpha2_;

        //- Mixture velocity
        const volVectorField& U_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- ...
        surfaceScalarField nHatf_;

        //- ...
        volScalarField K_;


    // Private Member Functions

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf
        );

        //- Re-calculate the interface curvature
        void calculateK();


public:

    // Constructors

        //- Construct from dictionary, volume fraction fields and mixture
        //  velocity
        interfaceProperties
        (
            const IOdictionary& dict,
            volScalarField& alpha1,
            volScalarField& alpha2,
            const volVectorField& U
        );

        //- Disallow default bitwise copy construction
        interfaceProperties(const interfaceProperties&) = delete;


    //- Destructor
    virtual ~interfaceProperties()
    {}


    // Member Functions

        inline const volVectorField& U() const
        {
            return U_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        tmp<volVectorField> n() const;

        tmp<volScalarField> sigmaK() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away from the interface.
        tmp<volScalarField> nearInterface() const;

        //- Interface fraction in a cell
        tmp<volScalarField::Internal> fraction() const;

        //- Correct the curvature
        void correct();

        //- Read phaseProperties dictionary
        bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interfaceProperties&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
