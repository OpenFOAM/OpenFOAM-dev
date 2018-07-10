/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::SingleKineticRateDevolatilisation

Description
    Single kinetic rate devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A1.exp(- E/R.T)

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
        R     = universal gas constant
        T     = temperature

    Usage:

        SingleKineticRateDevolatilisationCoeffs
        {
            volatileData
            (
                (CH4     12     0.5)   // (name A1 E)
                (CO2     12     0.5)   // (name A1 E)
            );

            volatileResidualCoeff 1e-6;
        }

\*---------------------------------------------------------------------------*/

#ifndef SingleKineticRateDevolatilisation_H
#define SingleKineticRateDevolatilisation_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class SingleKineticRateDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SingleKineticRateDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
        // Private data

            //- Specie name
            word name_;

            //- Activation energy
            scalar A1_;

            //- Pre-exponential factor
            scalar E_;


    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
                A1_(0.0),
                E_(0.0)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
                A1_(readScalar(is)),
                E_(readScalar(is))
            {}

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
                A1_(vd.A1_),
                E_(vd.E_)
            {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the activation energy
                scalar A1() const
                {
                    return A1_;
                }

                //- Return const access to the pre-exponential factor
                scalar E() const
                {
                    return E_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.A1_  >> vd.E_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.A1_ << token::SPACE
                    << vd.E_
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants

            //- List of volatile data - (name A1 E)
            List<volatileData> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;


public:

    //- Runtime type information
    TypeName("singleKineticRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        SingleKineticRateDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        SingleKineticRateDevolatilisation
        (
            const SingleKineticRateDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DevolatilisationModel<CloudType>>
            (
                new SingleKineticRateDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SingleKineticRateDevolatilisation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            label& canCombust,
            scalarField& dMassDV
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SingleKineticRateDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
