// Calculate viscous force: velocity-gradient * area * kinematic viscosity * density.
fishViscousForce	=	- U.boundaryField()[forcePatchID].snGrad()
				* mesh.magSf().boundaryField()[forcePatchID]
				* nu.value()
				* rhoF.value();

// Update the viscous force field.
viscousForce.boundaryField()[forcePatchID] == fishViscousForce;

// Calculate total force: viscous force + pressure * face area.
fishTotalForce 		=	  viscousForce
				+ p.boundaryField()[forcePatchID]
				* mesh.Sf().boundaryField()[forcePatchID];

// Update the total force field.
totalForce.boundaryField()[forcePatchID] == fishTotalForce;

// Calculate input power needed to impose undulation against the pressure field.
// Sum up the power exerted by each face of the patch. Power is pressure force & face velocity.
// Face velocity can be obtained directly from fluid velocity field because of non-slip condition.
// Since forces are directed outwards of the fluid domain because surface vectors point outwards face velocity has the be inverted.
// What about viscous forces?
inputPower = gSum(					(
			p.boundaryField()[forcePatchID] * mesh.Sf().boundaryField()[forcePatchID]	
			)
			& -U.boundaryField()[forcePatchID]		// Velocity works against fluid force.	
		 );
Info << "Undulation power = " << inputPower << endl;
		

// Does not work in parallel, probably gsum has to be used in any case...
// Calculate total power by global force on patch and swimming speed from inlet.
// Only drag forces have to be considered, i.e. the x-component of the forces has to be tested for positive or negative direction and only positive forces will be integrated.
vector force;
scalar dragForce = 0;
outputPower = 0;	// Works for one patch only!
forAll(fishViscousForce, forceI)
{	
	// Calculate force on current face.
	force = fishViscousForce[forceI] + p.boundaryField()[forcePatchID][forceI] * mesh.Sf().boundaryField()[forcePatchID][forceI];

	// Integrate positive forces.
	if ( force[0] > 0 )
		{
		outputPower += (force & vector(1,0,0));		//U.boundaryField()[inletPatchID]
		dragForce += force[0];
		}
}

// Does not work in parallel as well...
/*	dragCheck = fishTotalForce;
forAll(dragForce.boundaryField()[forcePatchID], forceI)
{
	if ( dragForce.boundaryField()[forcePatchID][forceI][0] <= 0 )
		{
		dragCheck[forceI] = vector(0,0,0);			
		Info << "check" << endl;
		}
	else
		{
		dragCheck[forceI] = vector(1,0,0);			
		}
}

Info << dragCheck << endl;

outputPower = gSum(
			dragForce.boundaryField()[forcePatchID] * dragCheck 
			) & vector(1,0,0);	// free stream velocity.*/
	

Info << "Output power = " << outputPower << endl;

// Calculate froude efficiency by output power / (output power + input power).
froudeEfficiency = outputPower / (outputPower + inputPower);
Info << "Froude efficiency = " << froudeEfficiency << endl;	

//*************************************** EoF *********************************************************
