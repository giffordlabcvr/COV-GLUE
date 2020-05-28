// Get variations
var master_ref = 'REF_MASTER_WUHAN_HU_1';
var master_alignment = 'AL_GISAID_CONSTRAINED';

//var codingfeatures = [ "ORF_1a", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6",
//                        "NSP7", "NSP8", "NSP9", "NSP10", "NSP11", "ORF_1ab", "NSP12",
//                        "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF_3a", "E", "M",
//                        "ORF_6", "ORF_7a", "ORF_7b", "ORF_8", "N", "ORF_10", "ORF_3b" ];

var codingfeatures = [  "ORF_6" ]; // DEV

// Get the reference sequence amino acid residues for each feature into map data structure
var refseqFeatureAaMap = initialise_refseq_feature_aa_map();

// Iterate through the coding features, retrieving variations that meet the criterion
_.each(codingfeatures, function(codingfeature) {

    glue.logInfo("Processing feature "+codingfeature+" in alignment "+master_alignment);

    // Get reference aminos for this feature from map
    var featureMap = refseqFeatureAaMap.codingfeature;
    
    // Use alignment to calculate aa frequencies in each listed coding feature
    var featureResultMap = {};

	glue.inMode("alignment/"+master_alignment, function(){
	
		var resultList = glue.tableToObjects(glue.command(["amino-acid", "frequency", "-r", master_ref, "-f", codingfeature]));		
		_.each(resultList,function(resultObj){
		
		    // Get differences from consensus
		    var codonLabel = resultObj.codon;
		    var aminoAcid  = resultObj.aminoAcid;
		    var numMembers = resultObj.numMembers;
		    var pctMembers = resultObj.pctMembers;
    		glue.logInfo("Processing position "+codonLabel+" - frequency of '"+aminoAcid+"' = "+pctMembers);
		    
		});

	});
	
	//  Identify sequences containing the variations we selected
	

});

// Get the reference sequence amino acid residues for listed coding features into a map data structure
function initialise_refseq_feature_aa_map() {

    var refseqFeatureResultMap = {}; 

	// Iterate through the coding features
	_.each(codingfeatures, function(codingfeature) {

		var labelledCodons;
		glue.inMode("reference/"+master_ref+"/feature-location/"+codingfeature, function(){
		
			var featureResultMap = {};
			var resultList = glue.tableToObjects(glue.command(["amino-acid"]));		
			_.each(resultList,function(resultObj){

				featureResultMap[resultObj.codonLabel] = resultObj;
				var codonLabel = resultObj.codonLabel;
				var aminoAcid  = resultObj.definiteAas;
				glue.logInfo("Setting AA at position '"+codonLabel+"' to be:'"+aminoAcid+"'");

			});

			refseqFeatureResultMap[codingfeature] = featureResultMap;
		});

	});

	return refseqFeatureResultMap;
}




