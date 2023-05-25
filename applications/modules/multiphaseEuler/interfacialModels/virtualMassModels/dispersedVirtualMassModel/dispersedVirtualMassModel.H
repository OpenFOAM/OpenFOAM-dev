/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::dispersedVirtualMassModel

Description
    Model for virtual mass between two phases where one phase can be considered
    dispersed in the other and the virtual mass therefore characterised by a
    virtual mass coefficient

SourceFiles
    dispersedVirtualMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedVirtualMassModel_H
#define dispersedVirtualMassModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "virtualMassModel.H"
#include "dispersedPhaseInterface.H"

namespace Foam
{
namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                    Class dispersedVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedVirtualMassModel
:
    public virtualMassModel
{
protected:

    // Protected data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    // Constructors

        // Construct from a dictionary and an interface
        dispersedVirtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~dispersedVirtualMassModel();


    // Member Functions

        //- Return the virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const = 0;

        //- Return the phase-intensive virtual mass coefficient Ki
        //  used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> Ki() const;

        //- Return the virtual mass coefficient K
        //  used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
