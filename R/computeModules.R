setClass("moduleWeb", representation(originalWeb="matrix", moduleWeb="matrix", orderA="vector", orderB="vector", modules="matrix"));

computeModules = function(web, deep = FALSE, deleteOriginalFiles=TRUE, steps=1000000) {

	result		= cM(web, depth=1, nrOfModule=1, ytop=1, xleft=1, ybottom=dim(web)[1], 
					  xright=dim(web)[2], prev_orderA=c(1:dim(web)[1]), prev_orderB=c(1:dim(web)[2]), 
					  modules=matrix(c(0, 1, c(1:sum(dim(web)))), 1), 
					  deepCompute=deep, delete=deleteOriginalFiles, steps=steps);
	result[[4]]	= result[[4]][order(result[[4]][,1]),];

	new("moduleWeb", originalWeb=web, moduleWeb=result[[1]], orderA=result[[2]], orderB=result[[3]], modules=result[[4]]);
}


# This function actually prepares the recursive computation of the modules and returns an object of class "moduleWeb"
cM = function(web, depth, nrOfModule, ytop, xleft, ybottom, xright, prev_orderA, prev_orderB, modules, deepCompute, delete, steps) {

	result = list();

	webName = paste("web-", depth, "-", nrOfModule, sep="");

	# write web to file
	web2edges(web[ytop:ybottom, xleft:xright], webName=webName);
	argv = c("identifyModules", "-filename", paste(webName, ".pairs", sep=""), "-bipartite", "-steps", as.integer(steps));
	argv = as.character(argv);
	argc = as.integer(length(argv));

	.C("identifyModules", argc, argv, PACKAGE="bipartite");
# because of unresolved issues, the dll needs to be unloaded/reloaded after each run.
	LIBS <- .dynLibs()
	bipLIB <- which(unlist(sapply(LIBS, function(x) x[1])) == "bipartite")
	IMpath <- 	LIBS[[bipLIB]][[2]] # absolute path on the system to the dll!!!
	dyn.unload(IMpath)
	dyn.load(IMpath)
	
	# read in data from result files
	data = readModuleData(webName, deleteOriginalFiles=delete);

	n_a		= dim(web)[1];
	n_b		= dim(web)[2];
	n		= n_a + n_b;

	orderAFile	= data[[1]];
	orderBFile	= data[[2]];
	modulesFile	= data[[3]];
	#infoFile	= data[[4]];

	# Permutation of the graph and therefore actualization of the permutation vectors is necessary
	# if more than one module are suggested
	if(nrow(modulesFile) > 1) {

		# Actualization of the permutation vectors
		tempA			= c(1:dim(web)[1]);
		tempB			= c(1:dim(web)[2]);
		tempA[ytop:ybottom]	= tempA[ytop:ybottom][orderAFile];
		tempB[xleft:xright]	= tempB[xleft:xright][orderBFile-length(orderAFile)];
		orderAFile		= prev_orderA[tempA];
		orderBFile		= prev_orderB[tempB];

		result[[1]] = web[tempA, tempB];
		result[[2]] = orderAFile;
		result[[3]] = orderBFile;

		# The matrix M containing the information about the identified modules is formatted in the following way:
		# Each row i contains the information about one certain module while
		# M[i,1] represents the nesting depth of the module,
		# M[i,2] = 1 iff the module is the last one to plot within its nesting module and 0 else,
		# and the M[i,j] = j-offset iff node j-offset is part of the module and 0 else.
		offset_M = 2;
		offset_modulesFile = 1;
		nrOfModules = dim(modulesFile)[1];

		M									= matrix(0, nrOfModules, (sum(dim(web))+offset_M));
		M[, 1]								= depth;
		M[nrOfModules, 2]						= 1;

		colvals = append(prev_orderA[ytop:ybottom], prev_orderB[xleft:xright]+n_a);
		rowvals = which(modulesFile[, 2:ncol(modulesFile)] > 0) %% nrow(modulesFile);
		rowvals[rowvals == 0] = nrow(modulesFile);
		if(length(colvals) == length(rowvals)) {
			for(i in 1:length(colvals)) {
				M[rowvals[i], colvals[i]+offset_M] = colvals[i];
			}
		}

		modulesFile							= M;
		result[[4]]							= rbind(modules, M);

		# Computation of potential modules nested within the ones found until now
		if(deepCompute) {
			order = append(orderAFile, (orderBFile+n_a));

			# Apply the recursive function cM(...) to each module
			for(i in 1:nrow(modulesFile)) {
				mod = modulesFile[i, (offset_M+1):ncol(modulesFile)];
				mod = mod[order];

				# Calculate the coordinates of the part of the web we are looking at
				j = n_a+1;
				while(mod[j] == 0) { j = j+1; }			# calculate x-coordinate of left lower corner of module border
				xleft_new = j - n_a;

				j = n_a;
				while(mod[j] == 0) { j = j-1; }			# calculate y-coordinate of left lower corner of module border
				ybottom_new = j;

				j = n;
				while(mod[j] == 0) { j = j-1; }			# calculate x-coordinate of right upper corner of module border
				xright_new = j - n_a;

				j = 1;
				while(mod[j] == 0) { j = j+1; }			# calculate y-coordinate of right upper corner of module border
				ytop_new = j;

				# An invocation of cM(...) is necessary only if there is the possibility to find more than one submodule the current module consists of
				if((ybottom_new - ytop_new)+1 > 1 && (xright_new - xleft_new)+1 > 1) {
					print(paste("Recursive invocation (depth: ", depth+1, ", module nr. ", i, ")", sep=""));
					result = cM(web=result[[1]], depth+1, nrOfModule=i, ytop=ytop_new, xleft=xleft_new, 
								ybottom=ybottom_new, xright=xright_new, prev_orderA=result[[2]], prev_orderB=result[[3]], 
								modules=result[[4]], deepCompute, delete, steps)
				}
			}
		}

	}
	else {
		result[[1]] = web;
		result[[2]] = prev_orderA;
		result[[3]] = prev_orderB;
		result[[4]] = modules;
	}

	result;

}


# This function reads in the data written down to the hard drive disk by the invoked C code and optionally triggers the deletion of the appropriate files
readModuleData = function(webName=NULL, deleteOriginalFiles=TRUE) {
	orderAFile = as.vector(as.matrix(read.table(paste(webName, ".ordA", sep=""), header=TRUE, sep="\t")));
	orderBFile = as.vector(as.matrix(read.table(paste(webName, ".ordB", sep=""), header=TRUE, sep="\t")));
	modulesFile = as.matrix(read.table(paste(webName, ".mod", sep=""), header=TRUE, sep="\t"));
	#infoFile = read.table(paste(webName, ".info", sep=""), header=TRUE, sep="\t");

	if(deleteOriginalFiles) deleteModuleData(webName);

	result = list();
	result[[1]] = orderAFile;
	result[[2]] = orderBFile;
	result[[3]] = modulesFile;
	#result[[4]] = infoFile;

	result;
}


# This function deletes the files written down to the hard drive disk by the invoked C code
deleteModuleData = function(webName=NULL) {
	unlink(paste(webName, ".pairs", sep=""));
	unlink(paste(webName, "-names.lut", sep=""));
	unlink(paste(webName, ".ordA", sep=""));
	unlink(paste(webName, ".ordB", sep=""));
	unlink(paste(webName, ".mod", sep=""));
	unlink(paste(webName, ".info", sep=""));
}


listModuleInformation = function(moduleWebObject) {

	if(isCorrectModuleWebObject(moduleWebObject)) {

		result	= list();

		web	= slot(moduleWebObject, "originalWeb");
		modules	= slot(moduleWebObject, "modules");

		n_a	= nrow(web);
		n_b	= ncol(web);

		offset	= 2;

		for(depth in unique(modules[,1])) {
			result[[depth+1]] = list();

			counter = 1;

			for(i in 1:nrow(modules)) {
				if(modules[i,1] == depth) {
					result[[depth+1]][[counter]]		= list();
					result[[depth+1]][[counter]][[1]]	= rownames(web)[modules[i,(offset+1):(n_a+offset)][modules[i,(offset+1):(n_a+offset)] > 0]];
					result[[depth+1]][[counter]][[2]]	= colnames(web)[(modules[i,(n_a+offset+1):(n_a+n_b+offset)][modules[i,(n_a+offset+1):(n_a+n_b+offset)] > 0])-n_a];

					counter = counter + 1;
				}
			}
		}

		result;
	}
}


printoutModuleInformation = function(moduleWebObject) {

	if(isCorrectModuleWebObject(moduleWebObject)) {

		modules	= slot(moduleWebObject, "modules");

		listOfModules = listModuleInformation(moduleWebObject);

		linebreak = "\n";

		a = linebreak;

		for(depth in unique(modules[,1])) {
			for(i in 1:length(listOfModules[[depth+1]])) {
				a = paste(a, "Depth: ", depth, linebreak, linebreak, "Nr of module: ", i, linebreak, linebreak, "Rownames: ", linebreak, sep="");
				for(j in 1:length(listOfModules[[depth+1]][[i]][[1]])) {
					a = paste(a, paste("\t", listOfModules[[depth+1]][[i]][[1]][j], sep=""), sep=linebreak);
				}
				a = paste(a, linebreak, linebreak, "Colnames: ", linebreak, sep="");
				for(j in 1:length(listOfModules[[depth+1]][[i]][[2]])) {
					a = paste(a, paste("\t", listOfModules[[depth+1]][[i]][[2]][j], sep=""), sep=linebreak);
				}
				a = paste(a, linebreak, linebreak, "__________________________", linebreak, linebreak, sep="");
			}
			a = paste(a, "__________________________", linebreak, linebreak, sep="");
		}

		cat(a);
	}
}