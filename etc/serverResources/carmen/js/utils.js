function createNavBarSkeleton(wrappingNavSelector){
	var rName = getRootName();
	addFixedTopNavSkeleton(wrappingNavSelector, "Carmen", "mainNav", "projectNav");
// 	addNavLink("#mainNav", "Home", "/" + rName, "#siteHomeLink");
//	addNavDrop("#mainNav", "Regions", "regionDrop");
//	addNavDrop("#mainNav", "Samples All Mips", "samplesDrop");
//	addNavDrop("#mainNav", "Pre-clustering Extraction Stats", "extractionDrop");
}

function populateNavBar(wrappingNavSelector, names){
	var rName = getRootName();
   	//get grouping names
//	var groupNames = names["mipRegions"];
//	var regionLinkPre = "/" + rName + "/showRegionInfo/";
//	var regionName = (names.hasOwnProperty("regionName")) ? names["regionName"] : "";
//	//get samples names
//	var sampName = (names.hasOwnProperty("sample")) ? names["sample"] : "";
//	var sampleNames = names["samples"];
//	var sampleLinkPre = "/" + rName + "/showOneSampAllMipData/";
//	var sampleExtractionLinkPre = "/" + rName + "/showInitialReadStatsPerSample/";
//	
//	
//    d3.select("#regionDrop")
//		.selectAll("li")
//		.data(groupNames)
//		.enter()
//			.append("li")
//				.attr("class", function(d){
//					if(d == regionName){
//						return "active";
//					}else{
//						return "";
//					}
//				})
//			.append("a")
//				.attr("href", function(d){ return regionLinkPre + d;})
//				.text(function(d){return d;});
//    
//	d3.select("#samplesDrop")
//		.selectAll("li")
//		.data(sampleNames)
//		.enter()
//			.append("li")
//				.attr("class", function(d){
//					if(d == sampName){
//						return "active";
//					}else{
//						return "";
//					}
//				})
//			.append("a")
//				.attr("href",function(d){return sampleLinkPre + d;} )
//				.text(function(d){return d;});
//	
//	d3.select("#extractionDrop")
//	.selectAll("li")
//	.data(sampleNames)
//	.enter()
//		.append("li")
//			.attr("class", function(d){
//				if(d == sampName){
//					return "active";
//				}else{
//					return "";
//				}
//			})
//		.append("a")
//			.attr("href",function(d){return sampleExtractionLinkPre + d;} )
//			.text(function(d){return d;});
   /**@todo add the target and region links again 
    * 			d3.select("#mipGeneLink").attr("href","/" + rName + "/showGeneInfo/" + geneName );
			d3.select("#mipTargetLink").attr("href","/" + rName + "/showMipInfo/" + mipName ).text(mipName);
    */
}


function createNavBar(wrappingNavSelector, names){
	createNavBarSkeleton(wrappingNavSelector);
	populateNavBar(wrappingNavSelector, names);
}

