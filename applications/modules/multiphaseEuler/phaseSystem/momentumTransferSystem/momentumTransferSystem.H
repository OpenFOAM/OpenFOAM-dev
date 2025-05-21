/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::momentumTransferSystem

Description
    Class which provides interfacial momentum transfer between a number of
    phases. Drag, virtual mass, lift, wall lubrication and turbulent dispersion
    are all modelled. Access is provided to a various forces and coefficients
    needed by the solver.

SourceFiles
    momentumTransferSystem.C

\*---------------------------------------------------------------------------*/

#ifndef momentumTransferSystem_H
#define momentumTransferSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class blendedDragModel;
class blendedVirtualMassModel;
class blendedLiftModel;
class blendedWallLubricationModel;
class blendedTurbulentDispersionModel;

/*---------------------------------------------------------------------------*\
                 Class momentumTransferSystem Declaration
\*---------------------------------------------------------------------------*/

class momentumTransferSystem
:
    public IOdictionary
{
    // Private Type Definitions

        typedef HashPtrTable
        <
            volScalarField,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > KdTable;

        typedef HashPtrTable
        <
            blendedDragModel,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > dragModelTable;

        typedef HashPtrTable
        <
            blendedVirtualMassModel,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > virtualMassModelTable;

        typedef HashPtrTable
        <
            blendedLiftModel,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > liftModelTable;

        typedef HashPtrTable
        <
            blendedWallLubricationModel,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > wallLubricationModelTable;

        typedef HashPtrTable
        <
            blendedTurbulentDispersionModel,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > turbulentDispersionModelTable;


    // Private Data

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- Cached drag coefficients for dragCorrs
        mutable KdTable Kds_;

        //- Drag models
        const dragModelTable dragModels_;

        //- Virtual mass models
        const virtualMassModelTable virtualMassModels_;

        //- Lift models
        const liftModelTable liftModels_;

        //- Wall lubrication models
        const wallLubricationModelTable wallLubricationModels_;

        //- Turbulent dispersion models
        const turbulentDispersionModelTable turbulentDispersionModels_;


    // Private Member Functions

        //- Create IO object for an optional constant/momentumTransfer file
        IOobject io(const phaseSystem&) const;

        //- Return the models dictionary
        template<class ModelType>
        const dictionary& modelsDict() const;

        //- Add a field to the given result. If the result is not valid then
        //  copy the field in as if the result was initially zero.
        void addTmpField
        (
            tmp<surfaceScalarField>& result,
            const tmp<surfaceScalarField>& field
        ) const;

        //- Invert the ADVs matrix inplace
        void invADVs(List<UPtrList<scalarField>>& ADVs) const;

        //- Invert the ADVs matrix inplace
        template<class GeoMesh>
        void invADVs
        (
            PtrList<PtrList<GeometricField<scalar, GeoMesh>>>& ADVs
        ) const;


public:

    //- Runtime type information
    TypeName("momentumTransferSystem");


    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Constructors

        //- Construct from a phase system
        momentumTransferSystem(const phaseSystem&);


    //- Destructor
    virtual ~momentumTransferSystem();


    // Member Functions

        //- Return the explicit force fluxes for the cell-based algorithm, that
        //  do not depend on phase mass/volume fluxes, and can therefore be
        //  evaluated outside the corrector loop. This includes things like
        //  lift, turbulent dispersion, and wall lubrication.
        PtrList<surfaceScalarField> Fs() const;

        //- As Fs, but for the face-based algorithm
        PtrList<surfaceScalarField> Ffs() const;

        //- Return the inverse of the central + drag + virtual mass
        //  coefficient matrix
        void invADVs
        (
            const PtrList<volScalarField>& As,
            PtrList<volVectorField>& HVms,
            PtrList<PtrList<volScalarField>>& invADVs,
            PtrList<PtrList<surfaceScalarField>>& invADVfs
        ) const;

        //- Return the inverse of the central + drag + virtual mass
        //  coefficient matrix
        PtrList<PtrList<surfaceScalarField>> invADVfs
        (
            const PtrList<surfaceScalarField>& Afs,
            PtrList<surfaceScalarField>& HVmfs
        ) const;

        //- Return the phase diffusivity
        //  divided by the momentum central coefficient
        virtual tmp<surfaceScalarField> alphaDByAf
        (
            const PtrList<volScalarField>& rAs
        ) const;

        //- Return the flux corrections for the cell-based algorithm. These
        //  depend on phase mass/volume fluxes, and must therefore be evaluated
        //  inside the corrector loop.
        PtrList<surfaceScalarField> ddtCorrs() const;

        //- Set the cell and faces drag correction fields
        void dragCorrs
        (
            PtrList<volVectorField>& dragCorrs,
            PtrList<surfaceScalarField>& dragCorrf
        ) const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
