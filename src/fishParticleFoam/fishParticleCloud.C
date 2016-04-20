/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "fishParticleCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(fishParticle, 0);
    defineTemplateTypeNameAndDebug(Cloud<fishParticle>, 0);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fishParticleCloud::fishParticleCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    bool readFields
)
:
    Cloud<fishParticle>(mesh, cloudName, false),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            "particleProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rhop_(dimensionedScalar(particleProperties_.lookup("rhop")).value()),
    e_(dimensionedScalar(particleProperties_.lookup("e")).value()),
    mu_(dimensionedScalar(particleProperties_.lookup("mu")).value())
{
    if (readFields)
    {
        fishParticle::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fishParticleCloud::move(const dimensionedVector& g, int check)
{

	// Particle Diameter
	scalar d = 1e-3;
if(check == 1){
	// Injector 1
	vector pos1 = vector(-0.5,0.1,0.05);
	//Set initial velocity vector
	vector vel1=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell1=mesh_.findCell(pos1);
	if(cell1>=0) {
	fishParticle* ptr= new fishParticle(*this,pos1,cell1,d,vel1);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 2
	vector pos2 = vector(-0.5,0.08,0.05);
	//Set initial velocity vector
	vector vel2=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell2=mesh_.findCell(pos2);
	if(cell2>=0) {
	fishParticle* ptr= new fishParticle(*this,pos2,cell2,d,vel2);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 3
	vector pos3 = vector(-0.5,0.06,0.05);
	//Set initial velocity vector
	vector vel3=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell3=mesh_.findCell(pos3);
	if(cell3>=0) {
	fishParticle* ptr= new fishParticle(*this,pos3,cell3,d,vel3);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 4
	vector pos4 = vector(-0.5,0.04,0.05);
	//Set initial velocity vector
	vector vel4=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell4=mesh_.findCell(pos4);
	if(cell4>=0) {
	fishParticle* ptr= new fishParticle(*this,pos4,cell4,d,vel4);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 5
	vector pos5 = vector(-0.5,0.02,0.05);
	//Set initial velocity vector
	vector vel5=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell5=mesh_.findCell(pos5);
	if(cell5>=0) {
	fishParticle* ptr= new fishParticle(*this,pos5,cell5,d,vel5);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 6
	vector pos6 = vector(-0.5,0,0.05);
	//Set initial velocity vector
	vector vel6=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell6=mesh_.findCell(pos6);
	if(cell6>=0) {
	fishParticle* ptr= new fishParticle(*this,pos6,cell6,d,vel6);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 7
	vector pos7 = vector(-0.5,-0.02,0.05);
	//Set initial velocity vector
	vector vel7=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell7=mesh_.findCell(pos7);
	if(cell7>=0) {
	fishParticle* ptr= new fishParticle(*this,pos7,cell7,d,vel7);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 7
	vector pos8 = vector(-0.5,-0.04,0.05);
	//Set initial velocity vector
	vector vel8=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell8=mesh_.findCell(pos8);
	if(cell8>=0) {
	fishParticle* ptr= new fishParticle(*this,pos8,cell8,d,vel8);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 8
	vector pos9 = vector(-0.5,-0.06,0.05);
	//Set initial velocity vector
	vector vel9=vector(1,0,0);

	// Find cell at specified injection position and add particle here
	label cell9=mesh_.findCell(pos9);
	if(cell9>=0) {
	fishParticle* ptr= new fishParticle(*this,pos9,cell9,d,vel9);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 9
	vector pos10 = vector(-0.5,-0.08,0.05);
	//Set initial velocity vector
	vector vel10=vector(1,0,0);
	// Find cell at specified injection position and add particle here
	label cell10=mesh_.findCell(pos10);
	if(cell10>=0) {
	fishParticle* ptr= new fishParticle(*this,pos10,cell10,d,vel10);
	Cloud<fishParticle>::addParticle(ptr);
	}

	// Injector 11
	vector pos11 = vector(-0.5,-0.1,0.05);
	//Set initial velocity vector
	vector vel11=vector(1,0,0);
	// Find cell at specified injection position and add particle here
	label cell11=mesh_.findCell(pos11);
	if(cell11>=0) {
	fishParticle* ptr= new fishParticle(*this,pos11,cell11,d,vel11);
	Cloud<fishParticle>::addParticle(ptr);
	}
}


    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> nuInterp(nu);

    fishParticle::trackData td(*this, rhoInterp, UInterp, nuInterp, g.value());

    Cloud<fishParticle>::move(td);
}


// ************************************************************************* //
