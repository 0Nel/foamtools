/*---------------------------------------------------------------------------------------------------*\
	File sampleData.H

	Description
	A file is created which contains calulated parameters as uTau, ReTau etc.

\*---------------------------------------------------------------------------------------------------*/


// Create a pointer to an output stream object.
autoPtr<OFstream> dataFilePtr;

// Create the file and fill the header.
if (dataFilePtr.empty())
{
	Info<< endl;
	Info<< "Creating Parameter file." << endl;

        // Check if this is the master application of a parallel run. Slaves are not allowed to write.
	if (Pstream::master())
	{

		// Define the name of the directory.
		word name = "dataSampling";

		// Get the start time of the current application run.
		word startTimeName = runTime.timeName(runTime.startTime().value());

		// Set up the directory to which the data will be written.
		fileName dataDir;

		// Check wether this is a parallel or a serial run and adjust file paths accordingly.
		if (Pstream::parRun())
		{
			// Put the data in undecomposed case (Note: gives problems for distributed data running)
			dataDir = runTime.path()/".."/name/startTimeName;
		}
		else
		{
			dataDir = runTime.path()/name/startTimeName;
		}

		// Create directory if does not exist.
		mkDir(dataDir);

		// Open new file at start up
		dataFilePtr.reset(new OFstream(dataDir/(name + ".dat")));

		// Add a header to the file.
		dataFilePtr()
			<< "#Time"	<< tab
			<< "dragForce"	<< tab 	
			<< "inputPower"	<< tab 	
			<< "outputPower"	<< tab 
			<< "froudeEfficiency"	<< endl;
        }
}

//*************************************** EoF *********************************************************
