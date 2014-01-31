function jsAvailClick()
{
	if (document.getElementById("avail-more-1").style.display == "none")
	{
		// change the title text
		document.getElementById("avail-link").innerHTML = "&#x25B2; Available data &#x25B2;"
		
		// display info
		document.getElementById("avail-more-1").style.display = "block";
		document.getElementById("avail-more-2").style.display = "block";
		document.getElementById("avail-more-3").style.display = "block";
	}
	else
	{
		// change the title text
		document.getElementById("avail-link").innerHTML = "&#x25BC; Available data &#x25BC;"
		
		// hide info
		document.getElementById("avail-more-1").style.display = "none";
		document.getElementById("avail-more-2").style.display = "none";
		document.getElementById("avail-more-3").style.display = "none";
	}
}
function jsEallClick()
{
	document.getElementById("iEmin").disabled = !document.getElementById("iEmin").disabled;
	document.getElementById("iEmax").disabled = !document.getElementById("iEmax").disabled;
	document.getElementById("idE").disabled   = !document.getElementById("idE").disabled;
}
function jsDataMOver(cur)
{
	for (var i = 1; i <= Nstates; i++)
	for (var j = 0; j <= Nstates; j++)
	{
		elem = document.getElementById("dat-" + i + "-" + j);
		
		if (elem.style.borderWidth == "4px")
			continue;
		else if (elem == cur)
			elem.style.borderWidth = "2px";
		else
			elem.style.borderWidth = "1px";
	}
}
function jsDataClick(cur)
{
	var li_list = document.getElementById("ali");
	var lf_list = document.getElementById("alf");
	var li = li_list.options[li_list.selectedIndex].value;
	var lf = lf_list.options[lf_list.selectedIndex].value;
	
	for (var i = 1; i <= Nstates; i++)
	for (var j = 0; j <= Nstates; j++)
	{
		elem = document.getElementById("dat-" + i + "-" + j);
		
		if (elem == cur)
		{
			elem.style.borderWidth = "4px";
			document.getElementById("datadescr").innerHTML = cur.getAttribute("data-data_" + li + "_" + lf);
		}
		else
		{
			elem.style.borderWidth = "1px";
		}
	}
}
function jsDataAngular()
{
	var li = document.getElementById("ali").selectedIndex;
	var lf = document.getElementById("alf").selectedIndex;
	
	var strli = document.getElementById("ali").options[li].value;
	var strlf = document.getElementById("alf").options[lf].value;
	
	jsDataClick(0);
	
	for (var i = 1; i <= Nstates; i++)
	{
		if (i <= li)
			document.getElementById("head-i-" + i).style.background = "#888888";
		else
			document.getElementById("head-i-" + i).style.background = "#880000";
		
		if (i <= lf)
			document.getElementById("head-f-" + i).style.background = "#888888";
		else
			document.getElementById("head-f-" + i).style.background = "#880000";
		
		for (var f = 0; f <= Nstates; f++)
		{
			var elem = document.getElementById("dat-" + i + "-" + f);
			
			if (i <= li || (f > 0 && f <= lf))
			{
				// set to "disabled" colour
				elem.style.background = "#888888";
			}
			else
			{
				// set to white (will work always)
				elem.style.background = "#FFFFFF";
				
				// set to the stored colour (only where defined; override white)
				var colour = elem.getAttribute("data-colour_"+ strli + "_" + strlf);
				elem.style.background = "rgb(" + (255 - colour) + "," + (255 - colour) + ",255)";
			}
		}
	}
}

function jsDataSet()
{
	var elem;
	
	// 	1s -> 1s
	elem = document.getElementById("dat-1-1");
	elem.dataset.colour_s_s = 160;
	elem.dataset.data_s_s = "1s &rarr; 1s<br/>(elastic scattering)<br/><br/>Energies [Ry]:<br/><ul><li>0.05&ndash;0.97</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 2s
	elem = document.getElementById("dat-1-2");
	elem.dataset.colour_s_s = 150;
	elem.dataset.data_s_s = "1s &rarr; 2s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.75&ndash;0.97</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 2p
	elem = document.getElementById("dat-1-2");
	elem.dataset.colour_s_p = 140;
	elem.dataset.data_s_p = "1s &rarr; 2p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.75&ndash;0.97</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 3s
	elem = document.getElementById("dat-1-3");
	elem.dataset.colour_s_s = 130;
	elem.dataset.data_s_s = "1s &rarr; 3s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.888&ndash;0.97</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 3p
	elem = document.getElementById("dat-1-3");
	elem.dataset.colour_s_p = 120;
	elem.dataset.data_s_p = "1s &rarr; 3p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.888&ndash;0.97</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 3d
	elem = document.getElementById("dat-1-3");
	elem.dataset.colour_s_d = 110;
	elem.dataset.data_s_d = "1s &rarr; 3d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.888&ndash;0.97</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 4s
	elem = document.getElementById("dat-1-4");
	elem.dataset.colour_s_s = 60;
	elem.dataset.data_s_s = "1s &rarr; 4s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.9375&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 4p
	elem = document.getElementById("dat-1-4");
	elem.dataset.colour_s_p = 60;
	elem.dataset.data_s_p = "1s &rarr; 4p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.9375&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 4d
	elem = document.getElementById("dat-1-4");
	elem.dataset.colour_s_d = 60;
	elem.dataset.data_s_d = "1s &rarr; 4d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.9375&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 4f
	elem = document.getElementById("dat-1-4");
	elem.dataset.colour_s_f = 60;
	elem.dataset.data_s_f = "1s &rarr; 4f<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.9375&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 5s
	elem = document.getElementById("dat-1-5");
	elem.dataset.colour_s_s = 30;
	elem.dataset.data_s_s = "1s &rarr; 5s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.96&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 5p
	elem = document.getElementById("dat-1-5");
	elem.dataset.colour_s_p = 30;
	elem.dataset.data_s_p = "1s &rarr; 5p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.96&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 5d
	elem = document.getElementById("dat-1-5");
	elem.dataset.colour_s_d = 30;
	elem.dataset.data_s_d = "1s &rarr; 5d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.96&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 5f
	elem = document.getElementById("dat-1-5");
	elem.dataset.colour_s_f = 30;
	elem.dataset.data_s_f = "1s &rarr; 5f<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.96&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 1s -> 5g
	elem = document.getElementById("dat-1-5");
	elem.dataset.colour_s_g = 30;
	elem.dataset.data_s_g = "1s &rarr; 5g<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.96&ndash;0.97 (needs larger grid)</li><li>1.2&ndash;40 (needs more partial waves)</li></ul>";
	
	// 2s -> 1s
	elem = document.getElementById("dat-2-1");
	elem.dataset.colour_s_s = 60;
	elem.dataset.data_s_s = "2s &rarr; 1s<br/>(de-excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.001&ndash;0.21</li></ul>";
	
	// 2s -> 2s
	elem = document.getElementById("dat-2-2");
	elem.dataset.colour_s_s = 60;
	elem.dataset.data_s_s = "2s &rarr; 2s<br/>(elastic scattering)<br/><br/>Energies [Ry]:<br/><ul><li>0.001&ndash;0.21</li></ul>";
	
	// 2s -> 2p
	elem = document.getElementById("dat-2-2");
	elem.dataset.colour_s_p = 60;
	elem.dataset.data_s_p = "2s &rarr; 2p<br/>(degenerated transition)<br/><br/>Energies [Ry]:<br/><ul><li>0.001&ndash;0.21</li></ul>";
	
	// 2s -> 3s
	elem = document.getElementById("dat-2-3");
	elem.dataset.colour_s_s = 60;
	elem.dataset.data_s_s = "2s &rarr; 3s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.138&ndash;0.21</li></ul>";
	
	// 2s -> 3p
	elem = document.getElementById("dat-2-3");
	elem.dataset.colour_s_p = 60;
	elem.dataset.data_s_p = "2s &rarr; 3p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.138&ndash;0.21</li></ul>";
	
	// 2s -> 3d
	elem = document.getElementById("dat-2-3");
	elem.dataset.colour_s_d = 60;
	elem.dataset.data_s_d = "2s &rarr; 3d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.138&ndash;0.21</li></ul>";
	
	// 2s -> 4s
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_s_s = 20;
	elem.dataset.data_s_s = "2s &rarr; 4s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2s -> 4p
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_s_p = 20;
	elem.dataset.data_s_p = "2s &rarr; 4p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2s -> 4d
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_s_d = 20;
	elem.dataset.data_s_d = "2s &rarr; 4d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2s -> 4f
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_s_f = 20;
	elem.dataset.data_s_f = "2s &rarr; 4f<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2p -> 1s
	elem = document.getElementById("dat-2-1");
	elem.dataset.colour_p_s = 60;
	elem.dataset.data_p_s = "2p &rarr; 1s<br/>(de-excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.001&ndash;0.21</li></ul>";
	
	// 2p -> 2p
	elem = document.getElementById("dat-2-2");
	elem.dataset.colour_p_p = 60;
	elem.dataset.data_p_p = "2p &rarr; 2p<br/>(elastic scattering)<br/><br/>Energies [Ry]:<br/><ul><li>0.001&ndash;0.21</li></ul>";
	
	// 2p -> 2s
	elem = document.getElementById("dat-2-2");
	elem.dataset.colour_p_s = 60;
	elem.dataset.data_p_s = "2p &rarr; 2s<br/>(degenerated transition)<br/><br/>Energies [Ry]:<br/><ul><li>0.001&ndash;0.21</li></ul>";
	
	// 2p -> 3s
	elem = document.getElementById("dat-2-3");
	elem.dataset.colour_p_s = 60;
	elem.dataset.data_p_s = "2p &rarr; 3s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.138&ndash;0.21</li></ul>";
	
	// 2p -> 3p
	elem = document.getElementById("dat-2-3");
	elem.dataset.colour_p_p = 60;
	elem.dataset.data_p_p = "2s &rarr; 3p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.138&ndash;0.21</li></ul>";
	
	// 2p -> 3d
	elem = document.getElementById("dat-2-3");
	elem.dataset.colour_p_d = 60;
	elem.dataset.data_p_d = "2p &rarr; 3d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.138&ndash;0.21</li></ul>";
	
	// 2p -> 4s
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_p_s = 20;
	elem.dataset.data_p_s = "2p &rarr; 4s<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2p -> 4p
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_p_p = 20;
	elem.dataset.data_p_p = "2p &rarr; 4p<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2p -> 4d
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_p_d = 20;
	elem.dataset.data_p_d = "2p &rarr; 4d<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
	
	// 2p -> 4f
	elem = document.getElementById("dat-2-4");
	elem.dataset.colour_p_f = 20;
	elem.dataset.data_p_f = "2p &rarr; 4f<br/>(excitation)<br/><br/>Energies [Ry]:<br/><ul><li>0.1875&ndash;0.21 (needs larger grid)</li></ul>";
}
