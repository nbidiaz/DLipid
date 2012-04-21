DLipid: Network inference and dynamic simulation application.

================================

This is a Java and Matlab application for lipid flux analysis and
dynamic simulations. The application uses time-course mass spectrometry 
data to automatically construct remodeling network and infer dynamic 
system using B-Spline algorithm.

DLipid is composed of two well-defined sub applications: 
FluxNetwork and DynamicSimulation.

The code is on http://nbidiaz.github.com/DLipid/.

=================================

This application can be used in linux, mac, and windows.
To use this tool, you must first download and unpack DLipid.

> tar -zxvf DLipid.tar.gz

> unzip DLipid.zip


You will also need to install Matlab:
http://www.mathworks.com/downloads/web_downloads.



# Step 1: Construct flux network and simulation
==================================
> cd bin
> java -jar DLipid.jar <inputFile> <minTh> [-m] [-n]
* Usage: java -jar DLipid.jar ../DATA/18_0-18_1.txt 0.3

inputFile:
	-- Notice: The default input folder is "../DATA/".
	   	   If only a filename is given, the application
	   	   interprets it as "../DATA/filename".
	   	   Use "./filename" for files in "bin" folder.
	-- Specification:
				time1		time2		...
		species1	mean	std	mean	std
		Species2	mean	std	mean	std
		...

minTh <FLOAT>:
	-- Significance score threshold for FluxNetwork application.
	-- Must be positive.

OPTIONS:	
	-m FLOAT    Maximum threshold. If this value is specified, multiple runs of
		    simulations will be performed using thresholds from
		    minTh to this number at step size 0.1.

		
	-n INT      Iteration number in B-Spline algorithm, default is 1000. The algorithm 
		    will stop at next-closest iteration where the error is min within last
		    four iterations. The algorithm will also stop when error change is
		    less than 0.0001 between subsequent steps.


Then you will see Matlab application starts after a few seconds. 
You can quit Matlab after calculations finish.


# Step 2: Retrieve flux network and simulation results.
======================================================

Go to output folder.

If only a filename is given in Step 1, default output folder is "../RESULTS/$inputFile_prefix/".
Otherwise an output folder "$inputFilePath/$inputFile_prefix/" will be created.

If you excite the example command at Step 1, you will see the following files in the output folder .

	"log_18_0-18_10.3.txt":  stores parameter values and information from dynamic simulations.
	
	"sim_18_0-18_10.3.jpg":  shows the fit of simulation to observed kinetics.

	"18_0-18_10.1.grapheme":   	 is a xml format file of inferred flux network.

	"18_0-18_10.3.sif":	 is a sif format file of reaction network.

	-- Specification:

		source1	parameter1	target1
		source2	parameter2	target2
		...	


# LipidGraphical: A graphical interface for flux analysis.
======================================================

This application allows vitalization of .grahpml network files.

java -jar FluxNetworkGUI.jar



# DynamicSimulation: 
======================================================

This application allows simulation with user-defined .sif reaction network file.

> cd bin
> matlab

* Usage: main_remodeling_sim('../DATA/18_0-18_1.txt','../RESULTS/18_0-18_1/18_0-18_10.3.sif');



AUTHOR:
Norberto Diaz Diaz <ndiaz@upo.es> and Lu Zhang <zhanglv@bc.edu>.
Copyright © 2012 Chuang Lab Boston College.