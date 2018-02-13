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

Class
    Foam::radiation::radiativeIntensityRay

Description
    Radiation intensity for a ray in a given direction

SourceFiles
    radiativeIntensityRay.C

\*---------------------------------------------------------------------------*/

#ifndef radiativeIntensityRay_H
#define radiativeIntensityRay_H

#include "absorptionEmissionModel.H"
#include "blackBodyEmission.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// Forward declaration of classes
class fvDOM;

/*---------------------------------------------------------------------------*\
                    Class radiativeIntensityRay Declaration
\*---------------------------------------------------------------------------*/

class radiativeIntensityRay
{
public:

    static const word intensityPrefix;


private:

    // Private data

        //- Refence to the owner fvDOM object
        const fvDOM& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Absorption/emission model
        const absorptionEmissionModel& absorptionEmission_;

        //- Black body
        const blackBodyEmission& blackBody_;

        //- Total radiative intensity / [W/m2]
        volScalarField I_;

        //- Total radiative heat flux on boundary
        volScalarField qr_;

        //- Incident radiative heat flux on boundary
        volScalarField qin_;

        //- Emitted radiative heat flux on boundary
        volScalarField qem_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Number of wavelengths/bands
        label nLambda_;

        //- List of pointers to radiative intensity fields for given wavelengths
        PtrList<volScalarField> ILambda_;

        //- Global ray id - incremented in constructor
        static label rayId;

        //- My ray Id
        label myRayId_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        radiativeIntensityRay(const radiativeIntensityRay&);

        //- Disallow default bitwise assignment
        void operator=(const radiativeIntensityRay&);


public:

    // Constructors

        //- Construct form components
        radiativeIntensityRay
        (
            const fvDOM& dom,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label lambda,
            const absorptionEmissionModel& absEmmModel_,
            const blackBodyEmission& blackBody,
            const label rayId
        );


    //- Destructor
    ~radiativeIntensityRay();


    // Member functions

        // Edit

            //- Update radiative intensity on i direction
            scalar correct();

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

            //- Add radiative intensities from all the bands
            void addIntensity();


        // Access

            //- Return intensity
            inline const volScalarField& I() const;

            //- Return const access to the boundary heat flux
            inline const volScalarField& qr() const;

            //- Return non-const access to the boundary heat flux
            inline volScalarField& qr();

            //- Return non-const access to the boundary incident heat flux
            inline volScalarField& qin();

            //- Return non-const access to the boundary emitted heat flux
            inline volScalarField& qem();

            //- Return const access to the boundary incident heat flux
            inline const volScalarField& qin() const;

            //- Return const access to the boundary emitted heat flux
            inline const volScalarField& qem() const;

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the number of bands
            inline scalar nLambda() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            //- Return the radiative intensity for a given wavelength
            inline const volScalarField& ILambda(const label lambdaI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radiativeIntensityRayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
