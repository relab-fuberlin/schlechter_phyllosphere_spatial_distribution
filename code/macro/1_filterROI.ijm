//	Coordinates of cells in images
	
// 	Set global measurements
run("Fresh Start");
run("Set Measurements...", "area min redirect=None decimal=3");
setForegroundColor(0, 0, 0);
setBackgroundColor(0, 0, 0);

//	Set the directory where the labeled mask images were obtained from omnipose and output directory
#@ File[] (label="Select channel directories", style="directories") listdirs
#@ File (label= "Choose OUTPUT/FILTER directory", style = "directory") out

for (i=0; i < listdirs.length; i++) {
	list = getFileList(listdirs[i]);
	for (j = 0; j < list.length; j++) {
		if (endsWith(list[j], ".png")){
		open(listdirs[i] + "/" + list[j]);
		id = getImageID();
		print("Analysing image " + list[j]);
		run("Set Scale...", "distance=11.0133 known=1 unit=micron");
		selectImage(id);
		run("Measure");
		maxPix = getResult("Max");			// Get the max pixel value
		x = Array.getSequence(maxPix+1);
		for (k = 1; k < x.length; k++) {
			setThreshold(x[k], x[k]);
			run("Analyze Particles...", "display clear add");
			val = getResult("Area");
			if (val > 5) {
				setColor(0);
				roiManager("Fill");
				}
			if (val >= 5) {
				setColor(x[k]);
				roiManager("Fill");	
				}
			roiManager("Deselect");
			roiManager("Delete");
		}
		name = getTitle();
		newname = replace(name, "_cp_masks", "");
		saveAs("png", out + "/filter-" + newname);
		}
		while (nImages()>0) {			// Close all images
          	selectImage(nImages());  
          	run("Close");
		} 
    }
}