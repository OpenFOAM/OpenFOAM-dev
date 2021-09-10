/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::TDACChemistryModel

Description
    Extends standardChemistryModel by adding the TDAC method.

    References:
    \verbatim
        Contino, F., Jeanmart, H., Lucchini, T., & D’Errico, G. (2011).
        Coupling of in situ adaptive tabulation and dynamic adaptive chemistry:
        An effective method for solving combustion in engine simulations.
        Proceedings of the Combustion Institute, 33(2), 3057-3064.

        Contino, F., Lucchini, T., D'Errico, G., Duynslaegher, C.,
        Dias, V., & Jeanmart, H. (2012).
        Simulations of advanced combustion modes using detailed chemistry
        combined with tabulation and mechanism reduction techniques.
        SAE International Journal of Engines,
        5(2012-01-0145), 185-196.

        Contino, F., Foucher, F., Dagaut, P., Lucchini, T., D’Errico, G., &
        Mounaïm-Rousselle, C. (2013).
        Experimental and numerical analysis of nitric oxide effect on the
        ignition of iso-octane in a single cylinder HCCI engine.
        Combustion and Flame, 160(8), 1476-1483.

        Contino, F., Masurier, J. B., Foucher, F., Lucchini, T., D’Errico, G., &
        Dagaut, P. (2014).
        CFD simulations using the TDAC method to model iso-octane combustion
        for a large range of ozone seeding and temperature conditions
        in a single cylinder HCCI engine.
        Fuel, 137, 179-184.
    \endverbatim

SourceFiles
    TDACChemistryModelI.H
    TDACChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef TDACChemistryModel_H
#define TDACChemistryModel_H

#include "standardChemistryModel.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"
#include "DynamicField.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class TDACChemistryModel
:
    public standardChemistryModel<ThermoType>
{
    // Private member data

        //- Switch to select performance logging
        Switch log_;

        //- Temporary simplified mechanism concentration field
        DynamicField<scalar> sc_;

        //- Temporary map from complete to simplified concentration fields
        //  c -> sc
        List<label> cTos_;

        //- Temporary map from simplified to complete concentration fields
        //  sc -> c
        DynamicList<label> sToc_;

        //- Mechanism reduction method
        autoPtr<chemistryReductionMethod<ThermoType>> mechRedPtr_;

        //- Mechanism reduction method reference
        chemistryReductionMethod<ThermoType>& mechRed_;

        bool mechRedActive_;

        //- Tabulation method
        autoPtr<chemistryTabulationMethod<ThermoType>> tabulationPtr_;

        //- Tabulation method reference
        chemistryTabulationMethod<ThermoType>& tabulation_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;


    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("TDAC");


    // Constructors

        //- Construct from thermo
        TDACChemistryModel(const fluidReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        TDACChemistryModel(const TDACChemistryModel&) = delete;


    //- Destructor
    virtual ~TDACChemistryModel();


    // Member Functions

        //- Create and return a TDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        //- Return a reference to the list of mass fraction fields
        inline const PtrList<volScalarField>& Y();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcdt
        ) const;


        // Chemistry model functions (overriding functions in
        // standardChemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);


        // ODE functions (overriding functions in standardChemistryModel to take
        // into account the variable number of species)

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& c,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;


        // Mechanism reduction access functions

            //- Return access to the mechanism reduction method
            //  Needed by ISAT
            inline const chemistryReductionMethod<ThermoType>& mechRed() const;

            //- Allow the reduction method to reset the number of species
            inline void setNSpecie(const label newNs);

            //- Return true if specie i is active
            inline bool active(const label i) const;

            //- Set specie i active
            inline void setActive(const label i);

            //- Return the index in the complete set of species
            //  corresponding to the index si in the simplified set of species
            inline label sToc(const label si) const;

            //- Return the index in the simplified set of species
            //  corresponding to the index si in the complete set of species
            inline label cTos(const label ci) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TDACChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TDACChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TDACChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
