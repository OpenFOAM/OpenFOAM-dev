/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "constSolidThermo.H"

/* * * * * * * * * * * * * * Protected Member Functions   * * * * * * * * * * */

template<class Type>
Foam::VolField<Type> Foam::constSolidThermo::readProperty
(
    const word& name,
    const dimensionSet& dimensions
) const
{
    const dictionary& propDict(subDict(name));
    const word propType(propDict.lookup("type"));

    if (propType == "uniform" || propType == "zonal")
    {
        VolField<Type> vtf
        (
            IOobject
            (
                phasePropertyName(name),
                mesh().time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensioned<Type>(dimensions, propDict.lookup<Type>("value"))
        );

        if (propType == "zonal")
        {
            const dictionary& zonesDict(propDict.subDict("zones"));

            Info<< "    Reading " << name << " for zones:" << endl;

            Field<Type>& vtfIf = vtf;

            forAllConstIter(dictionary, zonesDict, iter)
            {
                const word& zoneName = iter().keyword();
                const Type value(pTraits<Type>(iter().stream()));

                Info<< "        " << zoneName << " " << value << endl;

                const labelList& zoneCells = mesh().cellZones()[zoneName];

                // Set the internal field in the zone to the value
                forAll(zoneCells, i)
                {
                    vtfIf[zoneCells[i]] = value;
                }

                // Set the patch field in the zone to the value

                const fvBoundaryMesh& patches = mesh().boundary();
                const labelList& own = mesh().faceOwner();

                boolList cellInZone(mesh().nCells(), false);

                forAll(zoneCells, i)
                {
                    cellInZone[zoneCells[i]] = true;
                }

                forAll(patches, patchi)
                {
                    const fvPatch& pp = patches[patchi];

                    forAll(pp, patchFacei)
                    {
                        const label facei = pp.start() + patchFacei;

                        if (cellInZone[own[facei]])
                        {
                            vtf.boundaryFieldRef()[patchi][patchFacei] = value;
                        }
                    }
                }
            }

            Info << endl;
        }

        vtf.correctBoundaryConditions();

        return vtf;
    }
    else if (propType == "file")
    {
        return VolField<Type>
        (
            IOobject
            (
                phasePropertyName(name),
                mesh().time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are "
            << "'uniform', 'zonal' or 'file' for " << name
            << abort(FatalError);

        return VolField<Type>::null();
    }
}


// ************************************************************************* //
