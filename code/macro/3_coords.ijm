//	Coordinates of cells in images
// 	Set global measurements
run("Fresh Start");

//	Set the directory where the labeled mask images were obtained from omnipose and output directory
#@ File (label= "Choose FILTER directory", style = "directory") dir
#@ File (label= "Choose COORD directory", style = "directory") out

list = getFileList(dir); // Get a list of all the files that will be analysed

setBatchMode(true);
for (i = 0; i < list.length; i++){
	if (startsWith(list[i], "sub-")){	// Select PNG files from the list
		open(dir + "/" + list[i]);		// Open individual files sequentially
		id = getImageID();
		title = getTitle();
		newT = replace(title, ".png", "");
		newT = replace(newT, "sub-filter-", "");
		print(list[i]);
		run("Set Scale...", "distance=11.0133 known=1 unit=micron");
		run("Set Measurements...", "min redirect=None decimal=3");
		selectImage(id);
    	run("Measure");
		maxPix = getResult("Max");// Get the max pixel value
		x = Array.getSequence(maxPix+1);
		run("Clear Results");
		run("Set Measurements...", "area center display invert redirect=None decimal=3");
		for (j = 0; j < x.length; j++){	//Loops through each label
			selectImage(id);
			setThreshold(x[j], x[j]);
			run("Analyze Particles...", "size=0.05-Infinity display exclude");
		}
		for (k = 0; k < nResults; k++) {
			setResult("Label", k, newT + "_" + k);
				}
		updateResults();
		selectWindow("Results");
		saveAs("Results", out + "/" + newT + ".csv");
	}
}
setBatchMode(false);
run("Fresh Start");