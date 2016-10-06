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
function jsDataMOver()
{
	var cur = this;
	
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
function jsDataClick()
{
	var cur = this;
	
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
	
	jsDataClick();
	
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
