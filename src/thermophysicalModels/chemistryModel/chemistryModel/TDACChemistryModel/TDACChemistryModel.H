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

        label timeSteps_;

        // Mechanism reduction
        label NsDAC_;
        scalarField completeC_;
        scalarField simplifiedC_;
        Field<bool> reactionsDisabled_;
        Field<label> completeToSimplifiedIndex_;
        DynamicList<label> simplifiedToCompleteIndex_;
        autoPtr<chemistryReductionMethod<ThermoType>> mechRed_;

        // Tabulation
        autoPtr<chemistryTabulationMethod<ThermoType>> tabulation_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;

        //- Log file for the average time spent adding tabulated data
        autoPtr<OFstream> cpuAddFile_;

        //- Log file for the average time spent growing tabulated data
        autoPtr<OFstream> cpuGrowFile_;

        //- Log file for the average time spent retrieving tabulated data
        autoPtr<OFstream> cpuRetrieveFile_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        // Field containing information about tabulation:
        // 0 -> add (direct integration)
        // 1 -> grow
        // 2 -> retrieve
        volScalarField tabulationResults_;


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

        //- Return the number of chemistry evaluations (used by ISAT)
        inline label timeSteps() const;

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

            inline void setNsDAC(const label newNsDAC);

            inline void setNSpecie(const label newNs);

            inline scalarField& completeC();

            inline scalarField& simplifiedC();

            inline Field<bool>& reactionsDisabled();

            inline bool active(const label i) const;

            inline void setActive(const label i);

            inline DynamicList<label>& simplifiedToCompleteIndex();

            inline Field<label>& completeToSimplifiedIndex();

            inline const Field<label>& completeToSimplifiedIndex() const;

            inline autoPtr<chemistryReductionMethod<ThermoType>>& mechRed();

            inline tmp<volScalarField> tabulationResults() const;

            void setTabulationResultsAdd(const label celli);

            void setTabulationResultsGrow(const label celli);

            void setTabulationResultsRetrieve(const label celli);

            inline void resetTabulationResults();


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
