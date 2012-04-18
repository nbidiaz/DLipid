DLipid: Network inference and dynamic symulation application.

# Step 1:

start shell application
go to bin directory

> java -jar DLipid.jar inputFile minimumTH [maximumTH]

   1. inputFile: 
	- Path: if only the file name is written, default input folder  will be selected: "../DATA/inputFile"
	- Specificiaton:
				time1		time2		...
		species1	mean	std	mean	std
		Species2	mean	std	mean	std
		...

   2.  Thresold used by FluxNetworkApplciaton: goes from the minimum thresold to the maximum one with the step size 0.1.
	- maximumTH is optional. If maximum threshold is not indicated, only the minimum threshold will be computed
   * Example: java -jar DLipid.jar 14 0-14 0.txt 0.2

# Step 2:

go to "DestinationFolder"/Results/"inputFile", see output files
   
   1. "DestinationFolder" depend of the inputFile path. Using 14_0-14_0.txt as example
	- if it is inside of DATA folder, ./DLipid/Results/14_0-14_0/
	- ohterwise; "inputFilePath"/14_0-14_0/

   2. Out files:
	- .txt: LULU please describes 
	- .eps: LULU please describes
	- .paghml: this file is an xml file format which stored flux network
	- .sif: codification of the remodeling network. Specification:

	  source1	parameter1	target1
	  source2	parameter2	target2
	  ...
		