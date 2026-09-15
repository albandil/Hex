/*
 * Hex scattering database -- scripts of the web interface
 */

// colour of the "full data" cell (matches the accent colour of the stylesheet)
var accentRGB = [ 138, 17, 19 ];

function jsAvailClick()
{
	var hidden = (document.getElementById("avail-more-1").style.display == "none");

	// change the title text
	document.getElementById("avail-link").innerHTML = hidden
		? "&#x25B2; Available data &#x25B2;"
		: "&#x25BC; Available data &#x25BC;";

	// show/hide the info
	for (var i = 1; i <= 3; i++)
	{
		var elem = document.getElementById("avail-more-" + i);

		if (elem != null)
			elem.style.display = hidden ? "block" : "none";
	}
}

function jsEallClick()
{
	document.getElementById("iEmin").disabled = !document.getElementById("iEmin").disabled;
	document.getElementById("iEmax").disabled = !document.getElementById("iEmax").disabled;
	document.getElementById("idE").disabled   = !document.getElementById("idE").disabled;
}

// add or remove a class of an element
function jsSetClass(elem, name, on)
{
	var classes = elem.className.split(/\s+/);
	var out = [];

	for (var i = 0; i < classes.length; i++)
		if (classes[i] != "" && classes[i] != name)
			out.push(classes[i]);

	if (on)
		out.push(name);

	elem.className = out.join(" ");
}

// iterate over all data cells of the availability table
function jsForEachCell(callback)
{
	for (var i = 1; i <= Nstates; i++)
	for (var j = 0; j <= Nstates; j++)
	{
		var elem = document.getElementById("dat-" + i + "-" + j);

		if (elem != null)
			callback(elem);
	}
}

function jsDataMOver()
{
	var cur = this;

	jsForEachCell(function (elem) {
		jsSetClass(elem, "hov", elem == cur);
	});
}

function jsDataClick()
{
	var cur = this;

	var li_list = document.getElementById("ali");
	var lf_list = document.getElementById("alf");
	var li = li_list.options[li_list.selectedIndex].value;
	var lf = lf_list.options[lf_list.selectedIndex].value;

	jsForEachCell(function (elem) {
		jsSetClass(elem, "sel", elem == cur);
	});

	if (cur != null && cur.getAttribute)
	{
		var descr = cur.getAttribute("data-data_" + li + "_" + lf);

		document.getElementById("datadescr").innerHTML = (descr != null)
			? descr
			: "No data available for this transition.";
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
		jsSetClass(document.getElementById("head-i-" + i), "disabled", i <= li);
		jsSetClass(document.getElementById("head-f-" + i), "disabled", i <= lf);

		for (var f = 0; f <= Nstates; f++)
		{
			var elem = document.getElementById("dat-" + i + "-" + f);

			if (elem == null)
				continue;

			if (i <= li || (f > 0 && f <= lf))
			{
				// forbidden transition (angular momentum too large)
				jsSetClass(elem, "disabled", true);
				elem.style.background = "";
			}
			else
			{
				// shade the cell according to the amount of the available data
				var colour = Number(elem.getAttribute("data-colour_" + strli + "_" + strlf)) || 0;
				var rgb = [];

				for (var c = 0; c < 3; c++)
					rgb.push(Math.round(255 - colour * (255 - accentRGB[c]) / 255));

				jsSetClass(elem, "disabled", false);
				elem.style.background = "rgb(" + rgb.join(",") + ")";
			}
		}
	}
}
