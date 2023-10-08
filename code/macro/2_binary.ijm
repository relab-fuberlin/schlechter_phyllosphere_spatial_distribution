
// 	Set global measurements
run("Fresh Start");


//	Directories
#@ File (label="Choose OUTPUT/FILTER input directory", style="directory" ) input
#@ File (label= "Choose OUTPUT/BINARY output directory", style = "directory") output
#@ String (label = "channel", value = "filter-C0") C0
#@ String (label = "channel", value = "filter-C1") C1
//#@ String (label = "channel", value = "filter-C2") C2

dir = getFileList(input);

// This chunk creates a subset of files based on channels
n = dir.length;
dirC0 = newArray;
dirC1 = newArray;
//dirC2 = newArray;

for (i = 0; i < n; i++) {
	if (startsWith(dir[i], C0)){
		dirC0 = Array.concat(dirC0, dir[i]);
	}
	else if (startsWith(dir[i], C1)){
		dirC1 = Array.concat(dirC1, dir[i]);
	}
//	else if (startsWith(dir[i], C2)){
//		dirC2 = Array.concat(dirC2, dir[i]);
//	}
}

// Function
function makeBin() {
	getMinAndMax(min, max);
	setThreshold(1, max);
	run("Convert to Mask");
	run("Dilate");
	run("Dilate");
	run("Erode");
	run("Erode");
}

setBatchMode(true);

//This bit is to obtain subtracted (non-overlapping) images

for (i=0; i < dirC0.length; i++){	// SynComm2
// for (i=0; i < dirC0.length; i++){ // SynComm3
	
	// Channel C0
	open(input + "/" + dirC0[i]);
	imgC0 = getTitle();
	rootName = replace(imgC0, "filter-C0-","");
	saveAs("png", input + "/sub-" + imgC0);
	makeBin();
	imgC0 = replace(imgC0, "filter-C0-", "binary-C0-");
	saveAs("png", output + "/" + imgC0);

	// Channel C1
	open(input + "/" + dirC1[i]);
	imgC1 = getTitle();
	imageCalculator("Subtract", imgC1, imgC0);
	saveAs("png", input + "/sub-" + imgC1);
	makeBin();
	imgC1 = replace(imgC1, "filter-C1-", "binary-C1-");
	saveAs("png", output + "/" + imgC1);
	
	// Channel C2
//	open(input + "/" + dirC2[i]);
//	imgC2 = getTitle();
//	imageCalculator("OR", imgC1, imgC0);
//	imageCalculator("Subtract", imgC2, imgC1);
//	saveAs("png", input + "/sub-" + imgC2);
//	makeBin();
//	imgC2 = replace(imgC2, "filter-C2-", "binary-C2-");
//	saveAs("png", output + "/" + imgC2);
	
	// Total biomass (C0 + C1)
	imageCalculator("OR", imgC0, imgC1);
	makeBin();
	saveAs("png", output + "/" + "binary-total-" + rootName);
	
	// Total biomass (C0 + C1 + C2)
//	imageCalculator("OR", imgC2, imgC1);
//	makeBin();
//	saveAs("png", output + "/" + "binary-total-" + rootName);
	
	while (nImages()>0) {
          selectImage(nImages());  
          run("Close");
    }
}