/*---------------------------------------------------------------------------------------------------*\
	File writeSampleData.H

	Description
	Data is written to a file.

\*---------------------------------------------------------------------------------------------------*/

// In a parallel run only the master is allowed to write.
if (Pstream::master())
{
	dataFilePtr()
		<< runTime.value()
		<< tab
		<< dragForce
		<< tab	
		<< inputPower
		<< tab
		<< outputPower
		<< tab
		<< froudeEfficiency
		<< tab
		<< endl;
}

//*************************************** EoF *********************************************************
