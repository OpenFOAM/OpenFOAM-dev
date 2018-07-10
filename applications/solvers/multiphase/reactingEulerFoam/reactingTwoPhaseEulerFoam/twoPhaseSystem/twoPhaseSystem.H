/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::twoPhaseSystem

Description
    Class which solves the volume fraction equations for two phases.

SourceFiles
    twoPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseSystem_H
#define twoPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dragModel;
class virtualMassModel;

/*---------------------------------------------------------------------------*\
                         Class twoPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem
:
    public phaseSystem
{
private:

    // Private member functions

        //- Return the drag coefficient for phase pair
        virtual tmp<volScalarField> Kd(const phasePairKey& key) const = 0;

        //- Return the face drag coefficient for phase pair
        virtual tmp<surfaceScalarField> Kdf(const phasePairKey& key) const = 0;

        //- Return the virtual mass coefficient for phase pair
        virtual tmp<volScalarField> Vm(const phasePairKey& key) const = 0;


protected:

    // Protected data

        //- Flag to indicate that returned lists of fields are "complete"; i.e.,
        //  that an absence of force is returned as a constructed list of zeros,
        //  rather than a null pointer
        static const bool fillFields_ = true;

        //- Phase model 1
        phaseModel& phase1_;

        //- Phase model 2
        phaseModel& phase2_;


public:

    //- Runtime type information
    TypeName("twoPhaseSystem");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            twoPhaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        twoPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~twoPhaseSystem();


    // Selectors

        static autoPtr<twoPhaseSystem> New
        (
            const fvMesh& mesh
        );


    // Member Functions

        using phaseSystem::sigma;
        using phaseSystem::dmdts;

        //- Return phase model 1
        const phaseModel& phase1() const;

        //- Access phase model 1
        phaseModel& phase1();

        //- Return phase model 2
        const phaseModel& phase2() const;

        //- Access phase model 2
        phaseModel& phase2();

        //- Return the phase not given as an argument
        const phaseModel& otherPhase(const phaseModel& phase) const;

        //- Return the surface tension coefficient
        tmp<volScalarField> sigma() const;

        //- Return the drag coefficient
        tmp<volScalarField> Kd() const;

        //- Return the face drag coefficient
        tmp<surfaceScalarField> Kdf() const;

        //- Return the virtual mass coefficient
        tmp<volScalarField> Vm() const;

        //- Solve for the phase fractions
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
