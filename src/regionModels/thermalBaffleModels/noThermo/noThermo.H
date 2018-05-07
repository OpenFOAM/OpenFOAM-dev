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
    Foam::regionModels::thermalBaffleModels::noThermo

Description
    Dummy surface pyrolysis model for 'none'

SourceFiles
    noThermo.C

\*---------------------------------------------------------------------------*/

#ifndef noThermo_H
#define noThermo_H

#include "thermalBaffleModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermalBaffleModels
{

/*---------------------------------------------------------------------------*\
                          Class noThermo Declaration
\*---------------------------------------------------------------------------*/

class noThermo
:
    public thermalBaffleModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        noThermo(const noThermo&);

        //- Disallow default bitwise assignment
        void operator=(const noThermo&);


protected:

    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from type name and mesh
        noThermo(const word& modelType, const fvMesh& mesh);

        //- Construct from components and dict
        noThermo
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict

        );


    //- Destructor
    virtual ~noThermo();


    // Member Functions


         // Thermo properties

            //- Return const reference to the solidThermo
            virtual const solidThermo& thermo() const;


        // Fields

            //- Return the film specific heat capacity [J/kg/K]
            virtual const tmp<volScalarField> Cp() const;

            //- Return solid absortivity [1/m]
            virtual const volScalarField& kappaRad() const;

            //- Return the film mean temperature [K]
            virtual const volScalarField& T() const;

            //- Return density [Kg/m3]
            virtual const volScalarField& rho() const;

            //- Return thermal conductivity [W/m/K]
            virtual const volScalarField& kappa() const;


        // Evolution

            //- Pre-evolve film
            virtual void preEvolveRegion();

            //- Evolve the film equations
            virtual void evolveRegion();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalBaffleModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
