
$(document).ready(function() {

	var rName = getRootName();
	var gifLoading = prsentDivGifLoading();
	var width = 575 * 2,
    	height = 560 * 2;

	var widthSeAsia = 575 * 2,
		heightSeAsia = 560 * 2;

	addDiv("body", "topNav");
	createNavBar("#topNav", {});
	d3.select(".navbar-brand").style("font-size", "30px")

	var mainDiv = addDiv("body", "mainContent").style("margin", "20px");

	addH1 ("#mainContent", "Sequences").style("opacity", 0);
	addH1 ("#mainContent", "Sequences").style("margin-left", "10px");
	addDiv("#mainContent", "combinedSeqs");

	postJSON("/" + rName + "/getCombinedSeqs", {}).then(function (seqData) {
		var seqViewerFinal = new njhSeqView("#combinedSeqs", seqData);
		var sesUid = seqData["sessionUID"];
		setUpCloseSession(sesUid);
	}).catch(logRequestError).then(function(){
		//done loading sites
	});




//	var projection = d3.geo.armadillo()
//	    .scale(400 * 1.75)
//	    .translate([width / 2, height * .7])
//	    .parallel(20)
//	    .rotate([-10, 0])
//	    .center([45.47633,2.618265])
//	    .precision(.1);

	var projection = d3.geo.armadillo()
	  .scale(400 * 2.5)
	  .translate([width / 2, height * .7])
	  .parallel(20)
	  .rotate([-10, 0])
	  .center([5,-10])
	  .precision(.1);

	var projectionSeAsia = d3.geo.armadillo()
	  .scale(400 * 5)
	  .translate([widthSeAsia / 2, heightSeAsia * .7])
	  .parallel(20)
	  .rotate([-10, 0])
	  .center([85,5])
	  .precision(.1);

	var arc = d3.svg.arc()
		.innerRadius(0)
		.outerRadius(10);

	var pie = d3.layout.pie()
		.value(function(d) { return d.proportion; }); // pie function transform data into specific format


	var path = d3.geo.path()
	    .projection(projection);

	var graticule = d3.geo.graticule();
	addH1 ("#mainContent", "Global Diversity").style("margin-left", "10px");
	var topMapDiv = d3.select("#mainContent").append("div").style("margin", "10px").attr("id", "topMapDiv")
	var topPhylogramDiv = addDiv("#mainContent", "phylogramTop").style("margin", "10px").style("clear", "both");
	addH1 ("#phylogramTop", "Phylogram").style("margin-left", "10px");
	var topPhylogramInputNamesDiv = addDiv("#phylogramTop", "phylogramInputNames")
		.style("margin", "10px")
		//.attr("id", "phylogramInputNames")
		.style("float", "left")
		//.style("margin", "5px")
		.style("overflow", "scroll")
		.style("height", height + "px")
		.style("width", 250 + "px");
	addH1 ("#phylogramInputNames", "Pop Names");
	var topPhylogramInputNamesSvg = topPhylogramInputNamesDiv.append("svg");

	var topPhylogramCountryNamesDiv = addDiv("#phylogramTop", "phylogramCountryNames")
		.style("margin", "10px")
		//.attr("id", "phylogramCountryNames")
		.style("float", "left")
		//.style("margin", "5px")
		.style("overflow", "scroll")
		.style("height", 1200 + "px")
		.style("width", 250 + "px");
	addH1 ("#phylogramCountryNames", "Country Names");
	var topPhylogramCountryNamesSvg = topPhylogramCountryNamesDiv.append("svg");

	var topPhylogramRefNamesDiv = addDiv("#phylogramTop", "phylogramRefNames")
		.style("margin", "10px")
		//.attr("id", "phylogramRefNames")
		.style("float", "left")
		//.style("margin", "5px")
		.style("overflow", "scroll")
		.style("height", 1200 + "px")
		.style("width", 250 + "px");
	addH1 ("#phylogramRefNames", "Reference Names");
	var topPhylogramRefNamesSvg = topPhylogramRefNamesDiv.append("svg");

	var phylogramDiv = addDiv("#phylogramTop", "phylogram")
			.style("margin", "20px")
			.style("float", "left");


	var inputNamesDiv = topMapDiv.append("div")
			.attr("id", "inputNamesDiv")
			.style("float", "left")
			//.style("margin", "5px")
			.style("overflow", "scroll")
			.style("height", height + "px")
			.style("width", 400 + "px");

	var namesSvg = inputNamesDiv.append("svg");



	var mapDiv = topMapDiv.append("div").style("float", "left")
	var svg = mapDiv.append("svg")
	    .attr("width", width)
	    .attr("height", height);

	var defs = svg.append("defs");

	defs.append("path")
	    .datum({type: "Sphere"})
	    .attr("id", "sphere")
	    .attr("d", path);

	defs.append("clipPath")
	    .attr("id", "clip")
	  .append("use")
	    .attr("xlink:href", "#sphere");

	svg.append("use")
	    .attr("class", "stroke")
	    .attr("xlink:href", "#sphere");

	svg.append("use")
	    .attr("class", "fill")
	    .attr("xlink:href", "#sphere");

	svg.append("path")
	    .datum(graticule)
	    .attr("class", "graticule")
	    .attr("clip-path", "url(#clip)")
	    .attr("d", path);

	var pathSeAsia = d3.geo.path()
	  .projection(projectionSeAsia);

	var graticule = d3.geo.graticule();
	var mapDivSeAsia = topMapDiv.append("div").style("float", "left")
	var svgSeAsia = mapDivSeAsia.append("svg")
	  .attr("width", widthSeAsia)
	  .attr("height", heightSeAsia)
    .style("margin-left", "100px");

	var defsSeAsia = svgSeAsia.append("defs");

	defsSeAsia.append("path")
	  .datum({type: "Sphere"})
	  .attr("id", "sphere")
	  .attr("d", pathSeAsia);

	defsSeAsia.append("clipPath")
	  .attr("id", "clip")
	.append("use")
	  .attr("xlink:href", "#sphere");

	svgSeAsia.append("use")
	  .attr("class", "stroke")
	  .attr("xlink:href", "#sphere");

	svgSeAsia.append("use")
	  .attr("class", "fill")
	  .attr("xlink:href", "#sphere");

	svgSeAsia.append("path")
	  .datum(graticule)
	  .attr("class", "graticule")
	  .attr("clip-path", "url(#clip)")
	  .attr("d", pathSeAsia);

	getJSON("/" + rName + "/getWorldJson").then(function (world) {
		getJSON("/" + rName + "/getSitesJson").then(function (sites) {
				//console.log(sites);
			getJSON("/" + rName + "/getCombinedJson").then(function (combinedJson) {
				getJSON("/" + rName + "/getAllNames").then(function (allNames) {
					//console.log(combinedJson);
					//console.log(Object.keys(combinedJson));
					var combinedJsonInput = [];

					for(obj in combinedJson){
						//console.log(obj);
						//console.log(combinedJson[obj]);
						if(combinedJson[obj].fromInput){
							combinedJsonInput.push(combinedJson[obj])
						}
					}
					//console.log(combinedJsonInput);
					namesSvg
							.style("height", (combinedJsonInput.length * 30 + 20) + "px")
							.style("width", 400 + "px");
					var nameGroups = namesSvg
						.selectAll(".nameGroup")
						.data(combinedJsonInput)
						.enter()
							.append("g")
								.attr("class", "nameGroup")
								.on("mouseover", function(d){
									var currentId = d.originalName;

									topMapDiv.selectAll(".arc").transition().attr('opacity',function () {
											//console.log("hello");
											//console.log(d3.select(this).datum().data.name);
								      return (d3.select(this).datum().data.name === currentId) ? 1.0 : 0.1;
										}).attr('stroke',function () {
										//console.log("hello");
										//console.log(d3.select(this).datum().data.name);
								      return (d3.select(this).datum().data.name === currentId) ? "white" : "black";
									});

								})
								.on("mouseout", function(d){
									topMapDiv.selectAll(".arc").transition().attr('opacity',function () {
	//									console.log("hello");
	//									console.log(d3.select(this).datum().data.name);
	//						      return (d3.select(this).datum().data.name === currentId) ? 1.0 : 0.25;
										return 0.75
									}).attr('stroke',function () {
										//console.log("hello");
										//console.log(d3.select(this).datum().data.name);
								      return "black";
									});
							})
					nameGroups.append("circle")
								.attr("class", "nameCircle")
								.attr("stroke", "#000000")
								.attr("fill", function(d,i){ return d.color})
								.attr("cx", 15)
								.attr("cy", function(d,i){ return i * 30 + 10})
								.attr("r", 10);
					nameGroups.append("text")
						.attr("class", "nameText")
						.attr("stroke", function(d,i){
							if(!d.foundInPop){
								return "red";
							}else{
								return "black";
							}
						})
						.attr('font-family', 'Helvetica Neue, Helvetica, sans-serif')
						.attr("x", 15 + 15)
						.attr("y", function(d,i){ return i * 30 + 12})
						.text(function(d) {
							if(d.refSeq){
								return d.name + "(" + d.refName + ")";
							}else{
								return d.name;
							}
							 });

					topPhylogramInputNamesSvg
						.style("height", (combinedJsonInput.length * 30 + 20) + "px")
						.style("width", 400 + "px");
					var phyNameGroups = topPhylogramInputNamesSvg
					.selectAll(".phyNameGroup")
					.data(combinedJsonInput)
					.enter()
						.append("g")
							.attr("class", "phyNameGroup");
					phyNameGroups.append("circle")
								.attr("class", "phyNameCircle")
								.attr("stroke", "#000000")
								.attr("fill", function(d,i){
										return "black";
										//return d.color
									})
								.attr("cx", 15)
								.attr("cy", function(d,i){ return i * 30 + 10})
								.attr("r", 10);
					phyNameGroups.append("text")
						.attr("class", "phyNameText")
						.attr("stroke", function(d,i){
							if(!d.foundInPop){
								//return "red";
								return "black";
							}else{
								return "black";
							}
						})
						.attr('font-family', 'Helvetica Neue, Helvetica, sans-serif')
						.attr("x", 15 + 15)
						.attr("y", function(d,i){ return i * 30 + 12})
						.text(function(d) {
							if(d.refSeq){
								return d.name + "(" + d.refName + ")";
							}else{
								return d.name;
							}
							 });

					//country names
					topPhylogramCountryNamesSvg
						.style("height", (allNames["countries"].length * 30 + 20) + "px")
						.style("width", 400 + "px");
					var phyCountryNameGroups = topPhylogramCountryNamesSvg
						.selectAll(".phyCountryNameGroup")
						.data(allNames["countries"])
						.enter()
							.append("g")
								.attr("class", "phyCountryNameGroup");
					phyCountryNameGroups.append("circle")
								.attr("class", "phyCountryNameCircle")
								.attr("stroke", "#000000")
								.attr("fill", function(d,i){
									return allNames["countryColors"][d.country];
										//console.log(allNames["countryColors"][d.country])
										//return "#" + allNames["countryColors"][d.country]["hexStr_"]
										//return "black";
										//return d.color
									})
								.attr("cx", 15)
								.attr("cy", function(d,i){ return i * 30 + 10})
								.attr("r", 10);
					phyCountryNameGroups.append("text")
						.attr("class", "phyCountryNameText")
						.attr("stroke", function(d,i){
							return "black";
						})
						.attr('font-family', 'Helvetica Neue, Helvetica, sans-serif')
						.attr("x", 15 + 15)
						.attr("y", function(d,i){ return i * 30 + 12})
						.text(function(d) {
							//console.log(d);
							return d.country;
							 });
				//ref names
					topPhylogramRefNamesSvg
					.style("height", (allNames["refNames"].length * 30 + 20) + "px")
					.style("width", 400 + "px");

					var phyRefNameGroups = topPhylogramRefNamesSvg
							.selectAll(".phyRefNameGroup")
							.data(allNames["refNames"])
							.enter()
								.append("g")
									.attr("class", "phyRefNameGroup");
						phyRefNameGroups.append("circle")
									.attr("class", "phyRefNameCircle")
									.attr("stroke", "#000000")
									.attr("fill", function(d,i){
											return "#" + allNames["refColors"][d]["hexStr_"]
											//return allNames["refColors"][d]
											//return "black";
											//return d.color
										})
									.attr("cx", 15)
									.attr("cy", function(d,i){ return i * 30 + 10})
									.attr("r", 10);
						phyRefNameGroups.append("text")
							.attr("class", "phyRefNameText")
							.attr("stroke", function(d,i){
								return "black";
							})
							.attr('font-family', 'Helvetica Neue, Helvetica, sans-serif')
							.attr("x", 15 + 15)
							.attr("y", function(d,i){ return i * 30 + 12})
							.text(function(d) {
								//console.log(d);
								return d;
								 });

				  svg.insert("path", ".graticule")
				      .datum(topojson.feature(world, world.objects.land))
				      .attr("class", "land")
				      .attr("clip-path", "url(#clip)")
				      .attr("d", path);

				  svg.insert("path", ".graticule")
				      .datum(topojson.mesh(world, world.objects.countries, function(a, b) { return a !== b; }))
				      .attr("class", "boundary")
				      .attr("clip-path", "url(#clip)")
				      .attr("d", path);
				  svg.selectAll(".citiesGroups")
				  	.data(sites)
				  	.enter()
				  		.append("g")
					  	.attr("class", "citiesGroups")
					  	.attr("transform",function(d){ return "translate(" + projection( d.coords )[0] + ", " + projection( d.coords )[1] + ")";})
					  	.selectAll(".arc")
					  		.data(function(d){ return pie(d.proportions);})
					  		.enter()
					  			.append("g")
					  			.attr("class", "arc")
					  			.attr("opacity", ".75")
					  			.attr("stroke", "#000000")
					  			.append("path")
					  				.attr("d", function(d){
					  					var arcSite = d3.svg.arc()
										.innerRadius(0)
										.outerRadius(Math.sqrt(d3.select(d3.select(this.parentNode).node().parentNode).datum().scaledSize)/Math.PI)
					  					return arcSite(d);})
					  				.attr("fill", function(d){
					  					return d.data.color;});

				  svgSeAsia.insert("path", ".graticule")
			      .datum(topojson.feature(world, world.objects.land))
			      .attr("class", "land")
			      .attr("clip-path", "url(#clip)")
			      .attr("d", pathSeAsia);

				  svgSeAsia.insert("path", ".graticule")
				      .datum(topojson.mesh(world, world.objects.countries, function(a, b) { return a !== b; }))
				      .attr("class", "boundary")
				      .attr("clip-path", "url(#clip)")
				      .attr("d", pathSeAsia);
				  svgSeAsia.selectAll(".citiesGroups")
				  	.data(sites)
				  	.enter()
				  		.append("g")
					  	.attr("class", "citiesGroups")
					  	.attr("transform",function(d){ return "translate(" + projectionSeAsia( d.coords )[0] + ", " + projectionSeAsia( d.coords )[1] + ")";})
					  	.selectAll(".arc")
					  		.data(function(d){ return pie(d.proportions);})
					  		.enter()
					  			.append("g")
					  			.attr("class", "arc")
					  			.attr("stroke", "#000000")
					  			.attr("opacity", ".75")
					  			.append("path")
					  				.attr("d", function(d){
					  					//console.log(Math.sqrt(d3.select(d3.select(this.parentNode).node().parentNode).datum().scaledSize)/Math.PI);
					  					var arcSite = d3.svg.arc()
										.innerRadius(0)
										.outerRadius(Math.sqrt(d3.select(d3.select(this.parentNode).node().parentNode).datum().scaledSize)/Math.PI)
					  					return arcSite(d);})
					  				.attr("fill", function(d){
					  					return d.data.color;});
					  				;
		//			var currentId = "PfDd2";
		//			console.log(currentId);
		//			console.log(d3.selectAll(".arc"))
		//			d3.selectAll(".arc").transition().style('opacity',function () {
		//				console.log("hello");
		//				console.log(d3.select(this).datum().data.name);
		//	      return (d3.select(this).datum().data.name === currentId) ? 1.0 : 0.25;
		//			});

				  makeRequest({url:"/" + rName + "/getCombinedNwk", method:"GET"}).then(function (combinedNwk) {

				var newick = NewickParse(combinedNwk, combinedJson);
//		        var newickNodes = []
//		        function buildNewickNodes(node, callback) {
//		          newickNodes.push(node)
//		          if (node.branchset) {
//		            for (var i=0; i < node.branchset.length; i++) {
//		              buildNewickNodes(node.branchset[i])
//		            }
//		          }
//		        }
//		        buildNewickNodes(newick)
		        //console.log(newick);
				var phy = d3.phylogram.build('#phylogram', newick, {
		          width: 1500,
		          height: 4500,
		          useInputNames:true,
		          skipLabels:true
		        });

					}).catch(logRequestError).then(function(){
						//done loading sites
					});
				}).catch(logRequestError).then(function(){
					//done loading sites
				});
			}).catch(logRequestError).then(function(){
				//done loading sites
			});
		}).catch(logRequestError).then(function(){
			//done loading sites
		});
	}).catch(logRequestError).then(function(){
		//done loading
		removeAllDivGifLoading();
	});

});
