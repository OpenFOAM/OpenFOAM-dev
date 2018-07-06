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

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveMesh::printAllocated() const
{
    Pout<< "primitiveMesh allocated :" << endl;

    // Topology
    if (cellShapesPtr_)
    {
        Pout<< "    Cell shapes" << endl;
    }

    if (edgesPtr_)
    {
        Pout<< "    Edges" << endl;
    }

    if (ccPtr_)
    {
        Pout<< "    Cell-cells" << endl;
    }

    if (ecPtr_)
    {
        Pout<< "    Edge-cells" << endl;
    }

    if (pcPtr_)
    {
        Pout<< "    Point-cells" << endl;
    }

    if (cfPtr_)
    {
        Pout<< "    Cell-faces" << endl;
    }

    if (efPtr_)
    {
        Pout<< "    Edge-faces" << endl;
    }

    if (pfPtr_)
    {
        Pout<< "    Point-faces" << endl;
    }

    if (cePtr_)
    {
        Pout<< "    Cell-edges" << endl;
    }

    if (fePtr_)
    {
        Pout<< "    Face-edges" << endl;
    }

    if (pePtr_)
    {
        Pout<< "    Point-edges" << endl;
    }

    if (ppPtr_)
    {
        Pout<< "    Point-point" << endl;
    }

    if (cpPtr_)
    {
        Pout<< "    Cell-point" << endl;
    }

    // Geometry
    if (cellCentresPtr_)
    {
        Pout<< "    Cell-centres" << endl;
    }

    if (faceCentresPtr_)
    {
        Pout<< "    Face-centres" << endl;
    }

    if (cellVolumesPtr_)
    {
        Pout<< "    Cell-volumes" << endl;
    }

    if (faceAreasPtr_)
    {
        Pout<< "    Face-areas" << endl;
    }

}


void Foam::primitiveMesh::clearGeom()
{
    if (debug)
    {
        Pout<< "primitiveMesh::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    deleteDemandDrivenData(cellCentresPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(cellVolumesPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
}


void Foam::primitiveMesh::clearAddressing()
{
    if (debug)
    {
        Pout<< "primitiveMesh::clearAddressing() : "
            << "clearing topology"
            << endl;
    }

    deleteDemandDrivenData(cellShapesPtr_);

    clearOutEdges();

    deleteDemandDrivenData(ccPtr_);
    deleteDemandDrivenData(ecPtr_);
    deleteDemandDrivenData(pcPtr_);

    deleteDemandDrivenData(cfPtr_);
    deleteDemandDrivenData(efPtr_);
    deleteDemandDrivenData(pfPtr_);

    deleteDemandDrivenData(cePtr_);
    deleteDemandDrivenData(fePtr_);
    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(ppPtr_);
    deleteDemandDrivenData(cpPtr_);
}


void Foam::primitiveMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


// ************************************************************************* //
