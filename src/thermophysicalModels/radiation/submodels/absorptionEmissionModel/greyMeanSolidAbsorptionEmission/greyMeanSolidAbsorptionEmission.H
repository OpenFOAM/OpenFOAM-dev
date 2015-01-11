/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::radiation::greyMeanSolidAbsorptionEmission

Description
    greyMeanSolidAbsorptionEmission radiation absorption and emission
    coefficients for continuous phase

    The coefficients for the species in the Look up table have to be specified
    for use in moles x P [atm], i.e. (k[i] = species[i]*p*9.869231e-6).

    The coefficients for CO and soot or any other added are multiplied by the
    respective mass fraction being solved

    All the species in the dictionary need either to be in the look-up table or
    being solved. Conversely, all the species solved do not need to be included
    in the calculation of the absorption coefficient

    The names of the species in the absorption dictionary must match exactly the
    name in the look-up table or the name of the field being solved

SourceFiles
    greyMeanSolidAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef greyMeanSolidAbsorptionEmission_H
#define greyMeanSolidAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "solidThermo.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class greyMeanSolidAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class greyMeanSolidAbsorptionEmission
:
    public absorptionEmissionModel
{

private:

    // Private data

        //- Enumering of radiative properties
        enum radiativeProperties
        {
            absorptivity,
            emissivity
        };

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- SLG thermo package
        const solidThermo& thermo_;

        //- Hash table of species names
        HashTable<label> speciesNames_;

        //- Basic multicomponent mixture
        const basicMultiComponentMixture& mixture_;

        //- List of solid species data
        List<FixedList<scalar, 2> > solidData_;


    // Private member functions

       //- Calculate the volumetric fraction of Yj
       tmp<scalarField> X(const word specie) const;

       //- Calculate the property mixing
       tmp<volScalarField> calc(const label) const;


public:

    //- Runtime type information
    TypeName("greyMeanSolidAbsorptionEmission");


    // Constructors

        //- Construct from components
        greyMeanSolidAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~greyMeanSolidAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
