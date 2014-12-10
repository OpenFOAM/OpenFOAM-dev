/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::radiation::wideBandAbsorptionEmission

Description

    wideBandAbsorptionEmission radiation absorption and emission coefficients
    for continuous phase.

    All the bands should have the same number of species and have to be entered
    in the same order.

    There is no check of continuity of the bands. They should not ovelap or
    have gaps.

    The emission constant proportionality is specified per band (EhrrCoeff).

    The coefficients for the species in the lookup table have to be specified
    for use in moles x P [atm].i.e. (k[i] = species[i]*p*9.869231e-6).

    The coefficients for CO and soot or any other added are multiplied by the
    respective mass fraction being solved.

    The look Up table file should be in the constant directory.

    band dictionary:
    \verbatim
        band0
        {
            bandLimits (1.0e-6 2.63e-6);
            EhrrCoeff       0.0;
            species
            {
                CH4
                {
                    Tcommon         300.;
                    Tlow            300.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (.1 0 0 0 0 0);
                }
                CO2
                {
                    Tcommon         300.;
                    Tlow            300.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (.1 0 0 0 0 0);
                }
                H2O
                {
                    Tcommon         300.;
                    Tlow            300.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (.1 0 0 0 0 0);
                }
                Ysoot
                {
                    Tcommon         300.;
                    Tlow            300.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (.1 0 0 0 0 0);
                }
            }
        }
    \endverbatim


SourceFiles
    wideBandAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef wideBandAbsorptionEmission_H
#define wideBandAbsorptionEmission_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wideBandAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wideBandAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of species considered for absorptivity
        static const int nSpecies_ = 5;

        //- Maximum number of bands
        static const int maxBands_ = 10;

        //-  Absorption coefficients
        FixedList<FixedList<absorptionCoeffs, nSpecies_>, maxBands_> coeffs_;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Hash table with species names
        HashTable<label> speciesNames_;

        //- Indices of species in the look-up table
        FixedList<label, nSpecies_> specieIndex_;

        //- Bands
        FixedList<Vector2D<scalar>, maxBands_> iBands_;

        //- Proportion of the heat released rate emitted
        FixedList<scalar, maxBands_> iEhrrCoeffs_;

        //- Lookup table of species related to ft
        mutable interpolationLookUpTable<scalar> lookUpTable_;

        //- Thermo package
        const fluidThermo& thermo_;

        //- Bands
        label nBands_;

        //- Pointer list of species being solved involved in the absorption
        UPtrList<volScalarField> Yj_;

        // Total wave length covered by the bands
        scalar totalWaveLength_;


public:

    //- Runtime type information
    TypeName("wideBandAbsorptionEmission");


    // Constructors

        //- Construct from components
        wideBandAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wideBandAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }

        //- Lower and upper limit of band i
        inline const Vector2D<scalar>& bands(const label i) const
        {
            return iBands_[i];
        }

        //- Correct rays
        void correct
        (
            volScalarField& a_,
            PtrList<volScalarField>& aLambda
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
